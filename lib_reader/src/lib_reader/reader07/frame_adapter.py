# -*- coding:utf-8 -*-
"""Unified reader for the GRID 04/07/09 (hexprint text) payloads.

Decoding (L1) and processing (L2) are split:

* :func:`_read_sci_frames` / :func:`_read_tel_frames` flatten the parsed packets
  into one row per particle / per telemetry sample (faithful L1 frames);
* :func:`frames_to_processed` reproduces the legacy ``single_read_hex``
  ``(sci, tel)`` output from those frames.

The three payloads share one packet format; the only differences are two
per-payload constants (``internal_resistance`` and the ``iMon`` scale divisor).
"""
from pathlib import Path

import numpy as np

from ..packet_parser import parse_grid_data_new
from ..frame_io import load_hex_text
from ..parity_check import crc16_xmodem_nd
from ..l1_cache import get_l1_frames, get_l2_processed, apply_selection

_XML = str(Path(__file__).with_name("grid_packet.xml"))

_INTERNAL_FREQ = 24.05e6
_SUB = 43
_ROWS_PER_PACKET = _SUB + 1

# per-payload constants (only the science-relevant differences between 04/07/09)
_PARAMS = {
    "04": {"internal_resistance": 2.1, "imon_div": 2.0},
    "07": {"internal_resistance": 1.1, "imon_div": 1.0},
    "09": {"internal_resistance": 1.1, "imon_div": 1.0},
}


def _crc_check(buf, start, nbytes, crc_offset):
    """crc16-xmodem over ``buf[start:start+nbytes]`` vs 2 bytes at crc_offset."""
    cov = buf[start:start + nbytes].reshape(1, -1)
    calc = int(crc16_xmodem_nd(cov)[0])
    stored = int(buf[crc_offset]) * 256 + int(buf[crc_offset + 1])
    return calc == stored


def _empty_sci_frames():
    return {
        "packet_idx": np.array([], dtype=np.int64),
        "event_idx": np.array([], dtype=np.int8),
        "channel": np.array([], dtype=np.uint8),
        "uscount": np.array([], dtype=np.uint64),
        "amp": np.array([], dtype=np.uint16),
        "effective_count": np.array([], dtype=np.int64),
        "missing_count": np.array([], dtype=np.int64),
        "crc_check": np.array([], dtype=bool),
    }


def _empty_tel_frames():
    out = {
        "crc_check": np.array([], dtype=bool),
        "utc": np.array([], dtype=np.uint32),
        "uscount": np.array([], dtype=np.uint64),
    }
    for name in ("temp_sipm", "temp_adc", "v_mon", "i_mon"):
        for i in range(4):
            out[f"{name}{i}"] = np.array([], dtype=np.uint16)
    return out


def _read_sci_frames(buf):
    """Decode the science hex stream into one row per particle (primary + 43 sub)."""
    data, index = parse_grid_data_new(
        "", xml_file=_XML, data_tag="hex_sci_packet", endian="MSB",
        multi_evt=_SUB, multi_step=11, data=buf,
    )
    if index.size == 0 or np.asarray(data["channel"]).size == 0:
        return _empty_sci_frames()

    good = np.array([_crc_check(buf, int(s), 510, int(s) + 510) for s in index[:, 0]])
    n = int(np.asarray(data["channel"]).size)
    prim_ch = np.asarray(data["channel"]).reshape(n)
    prim_us = np.asarray(data["uscount"]).reshape(n)
    prim_amp = np.asarray(data["data_max"]).reshape(n)
    sub_ch = np.asarray(data["sub_channel"]).reshape(n, _SUB)
    sub_us = np.asarray(data["sub_uscount"]).reshape(n, _SUB)
    sub_amp = np.asarray(data["sub_amp"]).reshape(n, _SUB)

    channel = np.concatenate([prim_ch[:, None], sub_ch], axis=1).ravel()
    uscount = np.concatenate([prim_us[:, None], sub_us], axis=1).ravel()
    amp = np.concatenate([prim_amp[:, None], sub_amp], axis=1).ravel()

    return {
        "packet_idx": np.repeat(np.arange(n, dtype=np.int64), _ROWS_PER_PACKET),
        "event_idx": np.tile(np.arange(_ROWS_PER_PACKET, dtype=np.int8), n),
        "channel": channel,
        "uscount": uscount,
        "amp": amp,
        "effective_count": np.repeat(
            np.asarray(data["effective_count"], dtype=np.int64), _ROWS_PER_PACKET),
        "missing_count": np.repeat(
            np.asarray(data["missing_count"], dtype=np.int64), _ROWS_PER_PACKET),
        "crc_check": np.repeat(good.astype(bool), _ROWS_PER_PACKET),
    }


def _read_tel_frames(buf):
    """Decode the telemetry hex stream into one row per sample (7 per packet)."""
    data, index = parse_grid_data_new(
        "", xml_file=_XML, data_tag="hex_tel_packet", endian="MSB",
        multi_evt=7, multi_step=70, data=buf,
    )
    if index.size == 0 or np.asarray(data["utc"]).size == 0:
        return _empty_tel_frames()

    good = np.array([_crc_check(buf, int(s), 496, int(s) + 496) for s in index[:, 0]])
    out = {
        "utc": np.asarray(data["utc"], dtype=np.uint32),
        "uscount": np.asarray(data["uscount"], dtype=np.uint64),
        "crc_check": np.repeat(good.astype(bool), 7),
    }
    for name in ("temp_sipm", "temp_adc", "v_mon", "i_mon"):
        arr = np.asarray(data[name])
        for i in range(4):
            out[f"{name}{i}"] = arr[:, i]
    return out


def frames_to_processed(sci_frames, tel_frames, internal_resistance, imon_div):
    """Reproduce the legacy 04/07/09 ``(sci, tel)`` output from L1 frames."""
    freq = _INTERNAL_FREQ

    # ---- science: primary event + 43 sub-events, per channel ----------------
    # Drop packets that failed the per-frame CRC (the legacy reader kept only
    # the good packet index, so this filter is part of the scientific output).
    # The legacy loop then walks each packet in order (primary event first, then
    # its 43 sub-events), so per-channel order follows packet order; reproduce
    # that order exactly via a stable group-by on the channel index.
    keep = np.asarray(sci_frames["crc_check"], dtype=bool)
    all_ch = sci_frames["channel"][keep].astype(np.int64) + 1
    all_amp = sci_frames["amp"][keep]
    all_ts = sci_frames["uscount"][keep].astype(np.float64) / freq

    valid = (all_ch > 0) & (all_ch < 5)
    all_ch = all_ch[valid]
    all_amp = all_amp[valid].astype(np.int64)
    all_ts = all_ts[valid]

    order = np.argsort(all_ch, kind="stable")
    ch_sorted = all_ch[order]
    amp_sorted = all_amp[order]
    ts_sorted = all_ts[order]
    bounds = np.concatenate([np.searchsorted(ch_sorted, [1, 2, 3, 4]), [len(ch_sorted)]])
    amp = [amp_sorted[bounds[i]:bounds[i + 1]] for i in range(4)]
    uscount_evt = [ts_sorted[bounds[i]:bounds[i + 1]] for i in range(4)]

    # one value per packet was repeated on every particle row at L1
    crc_per_packet = np.asarray(sci_frames["crc_check"], dtype=bool)[::_ROWS_PER_PACKET]
    eff_per_packet = np.asarray(sci_frames["effective_count"])[::_ROWS_PER_PACKET]
    miss_per_packet = np.asarray(sci_frames["missing_count"])[::_ROWS_PER_PACKET]
    effective_count = np.asarray(eff_per_packet[crc_per_packet], dtype=np.int64)
    missing_count = np.asarray(miss_per_packet[crc_per_packet], dtype=np.int64)

    # ---- telemetry: 7 records per packet, 4 channels ------------------------
    # Same per-frame CRC filter as the science path.
    tkeep = np.asarray(tel_frames["crc_check"], dtype=bool)
    utc = np.asarray(tel_frames["utc"][tkeep], dtype=np.float64)
    uscount = np.asarray(tel_frames["uscount"][tkeep], dtype=np.float64) / freq

    temp_sipm = np.stack(
        [np.asarray(tel_frames[f"temp_sipm{i}"][tkeep], dtype=np.float64) for i in range(4)], axis=1)
    temp_adc = np.stack(
        [np.asarray(tel_frames[f"temp_adc{i}"][tkeep], dtype=np.float64) for i in range(4)], axis=1)
    v_mon = np.stack(
        [np.asarray(tel_frames[f"v_mon{i}"][tkeep], dtype=np.float64) for i in range(4)], axis=1)
    i_mon = np.stack(
        [np.asarray(tel_frames[f"i_mon{i}"][tkeep], dtype=np.float64) for i in range(4)], axis=1)

    temp_sipm = np.where(temp_sipm > 2048, (temp_sipm - 4096) / 16.0, temp_sipm / 16.0)
    temp_adc = np.where(temp_adc > 2048, (temp_adc - 4096) / 16.0, temp_adc / 16.0)
    v_mon = v_mon / 4096.0 * 3.3 * 11.0
    i_mon = i_mon / 4096.0 * 3.3 / imon_div
    bias = v_mon - i_mon * internal_resistance

    sci_extracted = {
        "amp": [np.asarray(amp[i]) for i in range(4)],
        "timestampEvt": [np.asarray(uscount_evt[i]) for i in range(4)],
        "timeCorrect": [],
        "sciNum": [np.asarray([]) for _ in range(4)],
        "effectiveCount": effective_count,
        "missingCount": missing_count,
    }
    tel_extracted = {
        "tempSipm": [temp_sipm[:, i] for i in range(4)],
        "tempAdc": [temp_adc[:, i] for i in range(4)],
        "vMon": [v_mon[:, i] for i in range(4)],
        "iMon": [i_mon[:, i] for i in range(4)],
        "bias": [bias[:, i] for i in range(4)],
        "timestamp": uscount,
        "utc": utc,
        "telNum": [],
    }
    return sci_extracted, tel_extracted


def _process(path, ver, internal_resistance, imon_div, select=None):
    buf = load_hex_text(path)
    frames = get_l1_frames(ver, ver, path, {
        "sci": lambda: _read_sci_frames(buf),
        "tl": lambda: _read_tel_frames(buf),
    })
    sci = frames["sci"]
    if select is not None:
        sci = apply_selection(ver, ver, path, {}, select[0], select[1], sci)
    return frames_to_processed(
        sci, frames["tl"], internal_resistance, imon_div)


def single_read_hex(path, internal_resistance, imon_div, ver="07", overwrite_cache=False, select=None):
    """Return ``(sci, tel)`` for a 04/07/09 hexprint file (legacy-shaped)."""
    params = {"internal_resistance": internal_resistance, "imon_div": imon_div}
    l2_params = dict(params)
    if select is not None:
        l2_params["select"] = select[0]
    return get_l2_processed(
        ver, ver, path, l2_params,
        lambda: _process(path, ver, internal_resistance, imon_div, select=select),
        overwrite=overwrite_cache,
    )


def single_read07(path, config=None, **kwargs):
    p = _PARAMS["07"]
    return single_read_hex(path, kwargs.get("internal_resistance", p["internal_resistance"]),
                           kwargs.get("imon_div", p["imon_div"]), ver="07",
                           overwrite_cache=kwargs.get("overwrite_cache", False),
                           select=kwargs.get("select"))


def single_read09(path, config=None, **kwargs):
    p = _PARAMS["09"]
    return single_read_hex(path, kwargs.get("internal_resistance", p["internal_resistance"]),
                           kwargs.get("imon_div", p["imon_div"]), ver="09",
                           overwrite_cache=kwargs.get("overwrite_cache", False),
                           select=kwargs.get("select"))


def single_read04(path, config=None, **kwargs):
    p = _PARAMS["04"]
    return single_read_hex(path, kwargs.get("internal_resistance", p["internal_resistance"]),
                           kwargs.get("imon_div", p["imon_div"]), ver="04",
                           overwrite_cache=kwargs.get("overwrite_cache", False),
                           select=kwargs.get("select"))
