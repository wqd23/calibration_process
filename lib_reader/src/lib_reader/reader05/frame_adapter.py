# -*- coding:utf-8 -*-
"""Unified reader for the GRID 03B / 05B payloads, split into L1 and L2.

* L1 decoders (:func:`_decode_sci_l1` / :func:`_decode_hk_l1` /
  :func:`_decode_tl_l1`) return the faithful raw fields of every science frame
  (one row per particle) and every HK/timeline sample.
* L2 (:func:`frames_to_processed`) reproduces the legacy ``(sci, tel)`` output:
  CRC filtering, run splitting, unit conversion, UTC fit and the per-file time
  cut.

The science packets are decoded by the shared ``packet_parser``
(``reader05/grid_packet.xml``); HK / timeline decoding lives in
:mod:`reader05.readout`.
"""
import functools
import re
from pathlib import Path

import numpy as np

from ..packet_parser import parse_grid_data_new
from ..frame_io import load_binary
from ..util import data_refactor
from ..l1_cache import get_l1_frames, get_l1_meta, get_l2_processed
from . import readout

_XML = str(Path(__file__).with_name("grid_packet.xml"))

_INTERNAL_FREQ = readout.INTERNAL_FREQ
_FEATURE_EVENTS = 20
_WF_SAMPLES = 256


# --------------------------------------------------------------------------- #
# L1 decoders (faithful raw frames)
# --------------------------------------------------------------------------- #
def _empty_sci(feature_mode):
    zi = np.zeros(0, dtype=np.int32)
    z32 = np.zeros(0, dtype=np.uint32)
    z64 = np.zeros(0, dtype=np.uint64)
    zb = np.zeros(0, dtype=bool)
    if feature_mode:
        return {
            "packet_idx": zi, "channel": z32, "event_number": z32,
            "sample_length": z32, "evt_timestamp": z64, "evt_data_sum": z64,
            "evt_data_max": z32, "evt_data_base": z32, "evt_crc": z32,
            "crc_check": zb,
        }
    return {
        "packet_idx": zi, "channel": z32, "event_number": z32, "sample_length": z32,
        "timestamp": z64, "data_sum": z64, "data_max": z32, "data_base": z32,
        "CRC": z32, "crc_check": zb,
        "waveform_data": np.zeros((0, _WF_SAMPLES), dtype=np.uint32),
    }


def _parse_sci_l1(sci_raw, feature_mode):
    buf = np.frombuffer(sci_raw, dtype=np.uint8)
    if feature_mode:
        d, _ = parse_grid_data_new(
            "", xml_file=_XML, data_tag="sci_ft_packet", endian="MSB",
            multi_evt=_FEATURE_EVENTS, multi_step=24, data=buf,
        )
        n = len(d["channel"])
        rep = functools.partial(np.repeat, repeats=_FEATURE_EVENTS)
        return {
            "packet_idx": np.repeat(np.arange(n, dtype=np.int32), _FEATURE_EVENTS),
            "channel": rep(d["channel"]),
            "event_number": rep(d["event_number"]),
            "sample_length": rep(d["sample_length"]),
            "evt_timestamp": d["evt_timestamp"],
            "evt_data_sum": d["evt_data_sum"],
            "evt_data_max": d["evt_data_max"],
            "evt_data_base": d["evt_data_base"],
            "evt_crc": d["CRC"],
            "crc_check": d["crc_check"],
        }
    d, _ = parse_grid_data_new(
        "", xml_file=_XML, data_tag="sci_wf_packet", endian="MSB", data=buf,
    )
    n = len(d["channel"])
    return {
        "packet_idx": np.arange(n, dtype=np.int32),
        "channel": d["channel"],
        "event_number": d["event_number"],
        "sample_length": d["sample_length"],
        "timestamp": d["timestamp"],
        "data_sum": d["data_sum"],
        "data_max": d["data_max"],
        "data_base": d["data_base"],
        "CRC": d["CRC"],
        "crc_check": d["crc_check"],
        "waveform_data": d["waveform_data"],
    }


def _naive_sci_bytes(raw_data, udp_pos):
    """Concatenate every UDP packet payload (no chunking, no resync reset).

    This is the faithful full byte stream: legacy chunked readout silently
    drops the frames that straddle a chunk boundary or a resync gap, and those
    drops are recorded in ``legacy_drop_frame_idx`` so L2 can reproduce them.
    """
    out = bytearray()
    for p in udp_pos:
        out.extend(raw_data[p + 11:p + 11 + readout.UDP_PACK_LEN])
    return bytes(out)


def _legacy_drop_packets(raw_data, udp_pos, feature_mode):
    """Packet indices (in the naive full stream) that legacy readout dropped.

    Replays the chunked extraction (``maxUdpReadout`` batches plus the >100
    packet-ID resync reset) and keeps the byte offsets of the frames it parsed;
    any naive frame whose start offset is absent is a dropped packet.
    """
    n = len(udp_pos)
    p = readout.UDP_PACK_LEN
    pkt_start = np.empty(n, dtype=np.int64)
    off = 0
    for i, pos in enumerate(udp_pos):
        pkt_start[i] = off
        off += len(raw_data[pos + 11:pos + 11 + p])

    naive = _naive_sci_bytes(raw_data, udp_pos)
    _, nidx = _parse_raw(naive, feature_mode)
    naive_pos = {int(s): i for i, (s, _e) in enumerate(nidx)}

    max_udp = readout.MAX_UDP_READOUT
    kept = set()
    first = 0
    ids = []
    bl = bytearray()
    for ipos in range(n):
        bl.extend(raw_data[udp_pos[ipos] + 11:udp_pos[ipos] + 11 + p])
        ids.append(raw_data[udp_pos[ipos] + 8] * 256 + raw_data[udp_pos[ipos] + 9])
        if len(ids) > 2 and (ids[-1] - ids[-2]) > 100:
            bl.clear()
            ids.clear()
            first = ipos + 1
        if (ipos + 1) % max_udp == 0 or ipos == n - 1:
            if bl:
                _, cidx = _parse_raw(bytes(bl), feature_mode)
                base = int(pkt_start[first]) if first < n else len(naive)
                for s, _e in cidx:
                    kept.add(base + int(s))
            bl = bytearray()
            ids = []
            first = ipos + 1
    return [i for s, i in naive_pos.items() if s not in kept]


def _parse_raw(buf, feature_mode):
    buf = np.frombuffer(buf, dtype=np.uint8)
    if feature_mode:
        return parse_grid_data_new(
            "", xml_file=_XML, data_tag="sci_ft_packet", endian="MSB",
            multi_evt=_FEATURE_EVENTS, multi_step=24, data=buf,
        )
    return parse_grid_data_new(
        "", xml_file=_XML, data_tag="sci_wf_packet", endian="MSB", data=buf,
    )


def _decode_sci_l1(path, feature_mode, no_udp):
    raw_data = load_binary(path).tobytes()
    if no_udp:
        return _parse_sci_l1(raw_data, feature_mode), []
    udp_pos = readout.findPackPos(raw_data, re.compile(readout.PATTERNS["udp"], re.S))
    frames = _parse_sci_l1(_naive_sci_bytes(raw_data, udp_pos), feature_mode)
    return frames, _legacy_drop_packets(raw_data, udp_pos, feature_mode)


def _decode_hk_l1(hk_file, ending):
    with open(hk_file, "rb") as fin:
        hk_raw = fin.read()
    extracter = {
        "x_ray": readout.extractHKData_raw,
        "normal": readout.extractHKData_normal_raw,
        "03b": readout.extractHKData_03b_raw,
    }
    d = extracter[ending](hk_raw)
    out = {"timestamp": np.asarray(d["timestamp"])}
    for name in ("bias", "iMon", "temp", "iSys"):
        arr = np.asarray(d[name])
        for i in range(4):
            out[f"{name}{i}"] = arr[i]
    return out


def _decode_tl_l1(timeline_file, ending):
    with open(timeline_file, "rb") as fin:
        tl_raw = fin.read()
    d = (readout.extractTimelineData_03b_raw(tl_raw) if ending == "03b"
         else readout.extractTimelineData_raw(tl_raw))
    return {"utc": d["utc"], "pps": d["pps"], "timestamp": d["timestamp"]}


# --------------------------------------------------------------------------- #
# L2 processing (reproduces the legacy (sci, tel))
# --------------------------------------------------------------------------- #
def _sci_process(sci_frames, feature_mode, drop_packets):
    split_run_time = readout.SPLIT_RUN_TIME
    freq = _INTERNAL_FREQ
    data_max = [[] for _ in range(4)]
    baseline = [[] for _ in range(4)]
    timestamp_evt = [[] for _ in range(4)]
    event_id = [[] for _ in range(4)]
    sci_num = [[] for _ in range(4)]
    empty_channel = []

    if len(drop_packets) > 0:
        drop = np.asarray(drop_packets, dtype=np.int64)
        sel = ~np.isin(np.asarray(sci_frames["packet_idx"]), drop)
    else:
        sel = np.ones(len(sci_frames["packet_idx"]), dtype=bool)
    d = {k: np.asarray(v)[sel] for k, v in sci_frames.items()}

    if feature_mode:
        n = len(d["crc_check"]) // _FEATURE_EVENTS
        ok = d["crc_check"].reshape(n, _FEATURE_EVENTS).all(axis=1)
        okexp = np.repeat(ok, _FEATURE_EVENTS)
        ts = d["evt_timestamp"][okexp].astype(np.float64) / freq
        channel = d["channel"][okexp].astype(np.uint8)
        ev_id = d["event_number"][okexp].astype(np.uint32)
        amplitude = d["evt_data_max"][okexp].astype(np.uint16)
        mean_baseline = d["evt_data_base"][okexp].astype(np.uint16)
    else:
        ok = d["crc_check"]
        ts = d["timestamp"][ok].astype(np.float64) / freq
        channel = d["channel"][ok].astype(np.uint8)
        ev_id = d["event_number"][ok].astype(np.uint32)
        amplitude = d["data_max"][ok].astype(np.uint16)
        mean_baseline = d["data_base"][ok].astype(np.uint16)

    scisection = 1
    q_sci = np.where(ts[:-1] > ts[1:] + split_run_time)[0]
    cur_sci_num = np.ones(len(ts)) * scisection
    if len(q_sci) > 0:
        last_sci_pos = 0
        for isci in range(len(q_sci)):
            scisection += 1
            cur_sci_num[last_sci_pos:q_sci[isci] + 1] = scisection
            last_sci_pos = q_sci[isci] + 1

    for ich in range(4):
        if len(np.where(channel == ich)[0]) == 0:
            if ich not in empty_channel:
                empty_channel.append(ich)
    for ich in range(4):
        if ich in empty_channel:
            continue
        q_ch = np.where(channel == ich)[0]
        data_max[ich].extend(list(amplitude[q_ch]))
        baseline[ich].extend(list(mean_baseline[q_ch]))
        timestamp_evt[ich].extend(list(ts[q_ch]))
        event_id[ich].extend(list(ev_id[q_ch]))
        sci_num[ich].extend(list(cur_sci_num[q_ch]))

    amp = [[] for _ in range(4)]
    for ich in range(4):
        if ich not in empty_channel:
            amp[ich] = np.array(data_max[ich]) - np.array(baseline[ich])
    return (
        np.array(amp, dtype=object),
        np.array(timestamp_evt, dtype=object),
        np.array(event_id, dtype=object),
        np.array(sci_num, dtype=object),
    )


def frames_to_processed(filename_no_path, sci_frames, hk_frames, tl_frames,
                        feature_mode, ending, drop_packets):
    split_run_time = readout.SPLIT_RUN_TIME

    amp, timestamp_evt, event_id, sci_num = _sci_process(sci_frames, feature_mode, drop_packets)

    # ---- HK ----------------------------------------------------------------
    raw_hk = {name: np.stack([np.asarray(hk_frames[f"{name}{i}"]) for i in range(4)])
              for name in ("bias", "iMon", "temp", "iSys")}
    raw_hk["timestamp"] = np.asarray(hk_frames["timestamp"])
    hk_data = readout._hk_convert(raw_hk)

    temp = hk_data["temp"]
    bias = hk_data["bias"]
    i_mon = hk_data["iMon"]
    timestamp = hk_data["timestamp"]
    i_sys = hk_data["iSys"]

    section = 1
    q_tel = np.where(np.array(timestamp)[:-1] > np.array(timestamp)[1:] + split_run_time)[0]
    tel_num = np.ones(len(timestamp)) * section
    if len(q_tel) > 0:
        last_tel_pos = 0
        for itel in range(len(q_tel)):
            section += 1
            tel_num[last_tel_pos:q_tel[itel] + 1] = section
            last_tel_pos = q_tel[itel] + 1

    temp = np.array(temp)
    i_mon = np.array(i_mon)
    bias = np.array(bias)
    timestamp = np.array(timestamp)
    i_sys = np.array(i_sys)
    tel_num = np.array(tel_num)

    # ---- timeline ----------------------------------------------------------
    tl_timestamp = np.asarray(tl_frames["timestamp"], dtype=np.float64)
    if ending != "03b":
        tl_timestamp = tl_timestamp * 100.
    utc = readout.getUTC(filename_no_path, tl_timestamp, np.asarray(tl_frames["utc"]),
                         timestamp, False)
    utc = np.array(utc)

    # ---- time cut (per-file, from cutFileRef) -------------------------------
    time_cut = readout.CUT_FILE_REF.get(filename_no_path)
    if time_cut is not None:
        if isinstance(time_cut, list):
            q1 = (timestamp >= 0) * (timestamp <= time_cut[1])
        else:
            q1 = (timestamp >= time_cut)
        timestamp = timestamp[q1]
        utc = utc[q1]
        temp = temp[:, q1]
        i_mon = i_mon[:, q1]
        bias = bias[:, q1]
        i_sys = i_sys[:, q1]
        tel_num = tel_num[q1]
        for ich in range(4):
            if isinstance(time_cut, list):
                q2 = (np.array(timestamp_evt[ich]) >= time_cut[0]) * (np.array(timestamp_evt[ich]) <= time_cut[1])
            else:
                q2 = np.array(timestamp_evt[ich]) >= time_cut
            timestamp_evt[ich] = np.array(timestamp_evt[ich])[q2]
            amp[ich] = np.array(amp[ich])[q2]
            event_id[ich] = np.array(event_id[ich])[q2]
            sci_num[ich] = np.array(sci_num[ich])[q2]

    sci_extracted = {
        "amp": amp,
        "timestampEvt": timestamp_evt,
        "eventID": event_id,
        "sciNum": sci_num,
    }
    tel_extracted = {
        "tempSipm": temp,
        "iMon": i_mon,
        "bias": bias,
        "iSys": i_sys,
        "timestamp": timestamp,
        "utc": utc,
        "telNum": tel_num,
    }
    return data_refactor(sci_extracted), data_refactor(tel_extracted)


def _sibling_paths(path):
    filename_no_path = re.split(r"\\|/", path)[-1]
    file_path = path.split(filename_no_path)[0]
    file_split = re.split(r"rundata|.dat", filename_no_path)
    hk_file = file_path + file_split[0] + "HK" + file_split[1] + ".dat"
    timeline_file = file_path + file_split[0] + "TimeLine" + file_split[1] + ".dat"
    return filename_no_path, hk_file, timeline_file


def _process(path, ver, reader, feature_mode, no_udp, ending):
    filename_no_path, hk_file, timeline_file = _sibling_paths(path)
    params = {"feature_mode": feature_mode, "no_udp": no_udp}
    sci_extra = {}

    def _decode_sci():
        frames, drop = _decode_sci_l1(path, feature_mode, no_udp)
        sci_extra["legacy_drop_frame_idx"] = drop
        return frames

    sci_frames = get_l1_frames(
        ver, reader, path, {"sci": _decode_sci}, params, extra_meta=sci_extra,
    )["sci"]
    drop_packets = get_l1_meta(ver, reader, path, params, "sci").get(
        "legacy_drop_frame_idx", [])
    hk_frames = get_l1_frames(
        ver, reader, hk_file, {"hk": lambda: _decode_hk_l1(hk_file, ending)},
        {"ending": ending},
    )["hk"]
    tl_frames = get_l1_frames(
        ver, reader, timeline_file, {"tl": lambda: _decode_tl_l1(timeline_file, ending)},
        {"ending": ending},
    )["tl"]
    return frames_to_processed(
        filename_no_path, sci_frames, hk_frames, tl_frames, feature_mode, ending,
        drop_packets)


def _read(path, ver, reader, feature_mode, no_udp, ending, overwrite_cache=False):
    params = {"feature_mode": feature_mode, "no_udp": no_udp, "ending": ending}
    return get_l2_processed(
        ver, reader, path, params,
        lambda: _process(path, ver, reader, feature_mode, no_udp, ending),
        overwrite=overwrite_cache,
    )


def single_read05b_normal(path, config=None, **kwargs):
    return _read(path, "05B", "normal", feature_mode=False, no_udp=True, ending="normal",
                 overwrite_cache=kwargs.get("overwrite_cache", False))


def single_read05b_xray(path, config=None, **kwargs):
    return _read(path, "05B", "xray", feature_mode=False, no_udp=True, ending="x_ray",
                 overwrite_cache=kwargs.get("overwrite_cache", False))


def single_read03b(path, config=None, **kwargs):
    return _read(path, "03B", "03b", feature_mode=True, no_udp=False, ending="03b",
                 overwrite_cache=kwargs.get("overwrite_cache", False))


def src_read03b(path, config=None, **kwargs):
    return _read(path, "03B", "03b-src", feature_mode=False, no_udp=False, ending="03b",
                 overwrite_cache=kwargs.get("overwrite_cache", False))
