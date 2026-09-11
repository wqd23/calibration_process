# -*- coding:utf-8 -*-
"""Unified frame-based reader for the GRID 04/07/09 (hexprint text) payloads.

These three payloads share one packet format and one reader; the only
differences are two per-payload constants (``internal_resistance`` and the
``iMon`` scale divisor).  The byte stream is decoded by the shared
``packet_parser`` (``reader07/grid_packet.xml``) and post-processed here into
the same ``(sci, tel)`` dicts the legacy ``gridBasicFunctions02.dataReadout``
produced, so the scientific output is unchanged.
"""
from pathlib import Path

import numpy as np

from ..packet_parser import parse_grid_data_new
from ..frame_io import load_hex_text
from ..parity_check import crc16_xmodem_nd

_XML = str(Path(__file__).with_name("grid_packet.xml"))

_INTERNAL_FREQ = 24.05e6

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


def _read_sci(path, buf):
    data, index = parse_grid_data_new(
        path, xml_file=_XML, data_tag="hex_sci_packet", endian="MSB",
        multi_evt=43, multi_step=11, data=buf,
    )
    if index.size == 0:
        return data, index
    good = np.array([
        _crc_check(buf, int(s), 510, int(s) + 510) for s in index[:, 0]
    ])
    return data, index[good]


def _read_tel(path, buf):
    data, index = parse_grid_data_new(
        path, xml_file=_XML, data_tag="hex_tel_packet", endian="MSB",
        multi_evt=7, multi_step=70, data=buf,
    )
    if index.size == 0:
        return data, index
    good = np.array([
        _crc_check(buf, int(s), 496, int(s) + 496) for s in index[:, 0]
    ])
    return data, index[good]


def single_read_hex(path, internal_resistance, imon_div):
    """Return ``(sci, tel)`` for a 04/07/09 hexprint file (legacy-shaped)."""
    buf = load_hex_text(path)

    sci, _ = _read_sci(path, buf)
    tel, _ = _read_tel(path, buf)

    freq = _INTERNAL_FREQ

    # ---- science: primary event + 43 sub-events, per channel ----------------
    # The legacy loop walks each packet in order (primary event first, then its
    # 43 sub-events), so per-channel order follows packet order.  Reproduce that
    # order exactly via a stable group-by on the channel index.
    n = sci["channel"].size
    sub = 43
    prim_ch = sci["channel"] + 1
    prim_amp = sci["data_max"]
    prim_ts = sci["uscount"] / freq
    sub_ch = sci["sub_channel"].reshape(n, sub) + 1
    sub_amp = sci["sub_amp"].reshape(n, sub)
    sub_ts = sci["sub_uscount"].reshape(n, sub) / freq

    all_ch = np.concatenate([prim_ch[:, None], sub_ch], axis=1).ravel()
    all_amp = np.concatenate([prim_amp[:, None], sub_amp], axis=1).ravel()
    all_ts = np.concatenate([prim_ts[:, None], sub_ts], axis=1).ravel()
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

    effective_count = np.asarray(sci["effective_count"], dtype=np.int64)
    missing_count = np.asarray(sci["missing_count"], dtype=np.int64)

    # ---- telemetry: 7 records per packet, 4 channels ------------------------
    utc = np.asarray(tel["utc"], dtype=np.float64)
    uscount = np.asarray(tel["uscount"], dtype=np.float64) / freq

    temp_sipm = np.asarray(tel["temp_sipm"], dtype=np.float64)
    temp_adc = np.asarray(tel["temp_adc"], dtype=np.float64)
    v_mon = np.asarray(tel["v_mon"], dtype=np.float64)
    i_mon = np.asarray(tel["i_mon"], dtype=np.float64)

    temp_sipm = np.where(temp_sipm > 2048, (temp_sipm - 4096) / 16.0, temp_sipm / 16.0)
    temp_adc = np.where(temp_adc > 2048, (temp_adc - 4096) / 16.0, temp_adc / 16.0)
    v_mon = v_mon / 4096.0 * 3.3 * 11.0
    i_mon = i_mon / 4096.0 * 3.3 / imon_div
    bias = v_mon - i_mon * internal_resistance

    # ---- legacy output shape (data_refactor applies 4-channel split) --------
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


def single_read07(path):
    p = _PARAMS["07"]
    return single_read_hex(path, p["internal_resistance"], p["imon_div"])


def single_read04(path):
    p = _PARAMS["04"]
    return single_read_hex(path, p["internal_resistance"], p["imon_div"])
