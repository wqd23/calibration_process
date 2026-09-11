# -*- coding:utf-8 -*-
"""Unified frame-based reader for the GRID 03B / 05B payloads.

The science packets (waveform and feature) are decoded by the shared
``packet_parser`` (``reader05/grid_packet.xml``); the HK / timeline extraction
and the UTC fit are the same modular functions the legacy ``dataReadout`` uses,
so the scientific output is unchanged.
"""
import re
from pathlib import Path

import numpy as np

from ..packet_parser import parse_grid_data_new
from ..frame_io import load_binary
from ..util import data_refactor
from . import readout

_XML = str(Path(__file__).with_name("grid_packet.xml"))

_INTERNAL_FREQ = readout.INTERNAL_FREQ


def _parse_sci_bytes(sci_raw, feature_mode):
    """Decode a science byte chunk into the flat event dict ``extractSciEvents``
    produced, in frame order with CRC filtering."""
    buf = np.frombuffer(sci_raw, dtype=np.uint8)

    if feature_mode:
        d, _ = parse_grid_data_new(
            "", xml_file=_XML, data_tag="sci_ft_packet", endian="MSB",
            multi_evt=20, multi_step=24, data=buf,
        )
        n = len(d["channel"])
        ok = d["crc_check"].reshape(n, 20).all(axis=1)
        okexp = np.repeat(ok, 20)
        return {
            "timestampEvt": d["evt_timestamp"][okexp].astype(np.float64) / _INTERNAL_FREQ,
            "channel": np.repeat(d["channel"][ok], 20).astype(np.uint8),
            "eventID": np.repeat(d["event_number"][ok], 20).astype(np.uint32),
            "amplitude": d["evt_data_max"][okexp].astype(np.uint16),
            "meanBaseline": d["evt_data_base"][okexp].astype(np.uint16),
        }

    d, _ = parse_grid_data_new(
        "", xml_file=_XML, data_tag="sci_wf_packet", endian="MSB", data=buf,
    )
    ok = d["crc_check"]
    return {
        "timestampEvt": d["timestamp"][ok].astype(np.float64) / _INTERNAL_FREQ,
        "channel": d["channel"][ok].astype(np.uint8),
        "eventID": d["event_number"][ok].astype(np.uint32),
        "amplitude": d["data_max"][ok].astype(np.uint16),
        "meanBaseline": d["data_base"][ok].astype(np.uint16),
    }


def _sci_process(path, feature_mode, no_udp):
    """Mirror ``dataReadout``'s chunked science loop: for UDP-wrapped payloads
    the stream is read in ``maxUdpReadout`` runs (reusing ``extractSciRawData``,
    so frames straddling a run boundary are lost exactly as legacy did)."""
    raw_data = load_binary(path).tobytes()

    if no_udp:
        chunks = [raw_data]
    else:
        udp_pos = readout.findPackPos(raw_data, re.compile(readout.PATTERNS["udp"], re.S))
        max_udp = readout.MAX_UDP_READOUT
        total_runs = int(float(len(udp_pos)) / float(max_udp)) + 1
        chunks = []
        last_pos = 0
        for _ in range(total_runs):
            chunks.append(readout.extractSciRawData(raw_data, udp_pos, max_udp, last_pos))
            last_pos += max_udp

    split_run_time = readout.SPLIT_RUN_TIME
    data_max = [[] for _ in range(4)]
    baseline = [[] for _ in range(4)]
    timestamp_evt = [[] for _ in range(4)]
    event_id = [[] for _ in range(4)]
    sci_num = [[] for _ in range(4)]
    empty_channel = []

    scisection = 1
    for sci_raw in chunks:
        sci_data = _parse_sci_bytes(sci_raw, feature_mode)
        ts = sci_data["timestampEvt"]
        q_sci = np.where(ts[:-1] > ts[1:] + split_run_time)[0]
        cur_sci_num = np.ones(len(ts)) * scisection
        if len(q_sci) > 0:
            last_sci_pos = 0
            for isci in range(len(q_sci)):
                scisection += 1
                cur_sci_num[last_sci_pos:q_sci[isci] + 1] = scisection
                last_sci_pos = q_sci[isci] + 1

        for ich in range(4):
            if len(np.where(sci_data["channel"] == ich)[0]) == 0:
                if ich not in empty_channel:
                    empty_channel.append(ich)
        for ich in range(4):
            if ich in empty_channel:
                continue
            q_ch = np.where(sci_data["channel"] == ich)[0]
            data_max[ich].extend(list(sci_data["amplitude"][q_ch]))
            baseline[ich].extend(list(sci_data["meanBaseline"][q_ch]))
            timestamp_evt[ich].extend(list(sci_data["timestampEvt"][q_ch]))
            event_id[ich].extend(list(sci_data["eventID"][q_ch]))
            sci_num[ich].extend(list(cur_sci_num[q_ch]))

    amp = [[] for _ in range(4)]
    for ich in range(4):
        if ich not in empty_channel:
            amp[ich] = np.array(data_max[ich]) - np.array(baseline[ich])
    amp = np.array(amp, dtype=object)
    timestamp_evt = np.array(timestamp_evt, dtype=object)
    event_id = np.array(event_id, dtype=object)
    sci_num = np.array(sci_num, dtype=object)
    return amp, timestamp_evt, event_id, sci_num


def _data_readout(path, feature_mode, no_udp, ending):
    filename_no_path = re.split(r"\\|/", path)[-1]
    file_path = path.split(filename_no_path)[0]
    file_split = re.split(r"rundata|.dat", filename_no_path)
    hk_file = file_path + file_split[0] + "HK" + file_split[1] + ".dat"
    timeline_file = file_path + file_split[0] + "TimeLine" + file_split[1] + ".dat"

    amp, timestamp_evt, event_id, sci_num = _sci_process(path, feature_mode, no_udp)
    split_run_time = readout.SPLIT_RUN_TIME

    # ---- HK ----------------------------------------------------------------
    with open(hk_file, "rb") as fin:
        hk_raw = fin.read()
    hk_extracter = {
        "x_ray": readout.extractHKData,
        "normal": readout.extractHKData_normal,
        "03b": readout.extractHKData_03b,
    }
    hk_data = hk_extracter[ending](hk_raw)

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
    with open(timeline_file, "rb") as fin:
        tl_raw = fin.read()
    if ending == "03b":
        tl_data = readout.extractTimelineData_03b(tl_raw)
    else:
        tl_data = readout.extractTimelineData(tl_raw)

    utc = readout.getUTC(filename_no_path, tl_data["timestamp"], tl_data["utc"], timestamp, False)
    utc = np.array(utc)

    # ---- time cut (per-file, from cutFileRef) -------------------------------
    time_cut = None
    if filename_no_path in readout.CUT_FILE_REF:
        time_cut = readout.CUT_FILE_REF[filename_no_path]
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


def single_read05b_normal(path):
    return _data_readout(path, feature_mode=False, no_udp=True, ending="normal")


def single_read05b_xray(path, config):
    return _data_readout(path, feature_mode=False, no_udp=True, ending="x_ray")


def single_read03b(path, config):
    return _data_readout(path, feature_mode=True, no_udp=False, ending="03b")


def src_read03b(path, config):
    return _data_readout(path, feature_mode=False, no_udp=False, ending="03b")
