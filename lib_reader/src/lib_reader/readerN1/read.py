# -*- coding:utf-8 -*-
"""
Reader for GRID N1 payload data.

Differences from reader12 (GRID 12B):
- N1 .event.dat files contain 584-byte grid1x_ft_packet (12B: 528-byte);
  38 events per packet with a 14-byte per-event stride (extra data_ccm field).
  The CLYC dataset instead records 1080-byte grid1x_wf_packet (512-point
  waveform variant, grid_packet_512wave.xml) -- the onboard-computed
  data_max/data_base/data_sum fields are still present, so amp is computed
  the same way; pass mode="wf" for CLYC files (chunked parse + CRC filter)
- N1 .hk files are the 178-byte `hk_packet` defined in yingtian_packet.xml
  (same 0x1a2b3c4d head magic as the 187-byte grid1x_hk_packet of 12B, but a
  different field layout)
- unlike 12B ground data, N1 sci packets carry a REAL utc (seconds), so sci
  events are matched to HK records by nearest utc and HK epochs are cut by
  time, not by (bias, temperature) binning
- file naming is `{temp}-{bias}-{idx}.event.dat` for single-bias points
  (e.g. m20C-265-125.event.dat) and `{temp}-{lo}-{hi}-{idx}.event.dat` for
  bias-scan files that cover several bias points in one run
  (e.g. 0C-265-290-139.event.dat); HK pairing is by the leading
  `{temp}-{bias...}` stem prefix to `{prefix}-ecu_{NNN}.hk`
- scan files are split into per-bias segments by UTC gaps (same approach as
  the gridN_cali L0 notebooks): cut at the midpoints of gaps with no events
  within +-2 s; pass `seg_bias` (bias code, e.g. 285) to select the segment
  whose mean bias matches within 0.1 V
"""
from pathlib import Path

import numpy as np
from addict import Dict
from cachier import cachier

from .parse_grid_data import parse_grid_data_new
from .parse_grid_iter import parse_grid_data_iter

VOL_SET = [265, 270, 275, 280, 283, 285, 287, 290]


@cachier(cache_dir=Path(".cache") / "N1", separate_files=True)
def readSci(path, mode="ft"):
    if mode == "ft":
        return Dict(
            parse_grid_data_new(
                path, data_tag="grid1x_ft_packet", endian="MSB", multi_evt=38, multi_step=14
            )[0]
        )
    elif mode == "wf":
        # CLYC waveform runs: 1080-byte packets, chunked parse, drop CRC fails
        parts = []
        for chunk in parse_grid_data_iter(
            path,
            data_tag="grid1x_wf_packet",
            endian="MSB",
            xml_file=str(Path(__file__).parent / "grid_packet_512wave.xml"),
            packet_len=1080,
            chunk_events=50000,
        ):
            d = chunk["data"]
            ok = np.asarray(d["crc_check"], dtype=bool)
            parts.append({k: np.asarray(v)[ok] for k, v in d.items()
                          if np.asarray(v).shape[:1] == ok.shape[:1]})
        merged = {}
        for k in parts[0]:
            merged[k] = np.concatenate([p[k] for p in parts])
        return Dict(merged)
    else:
        raise ValueError("Invalid mode. Use 'ft' or 'wf'.")


@cachier(cache_dir=Path(".cache") / "N1", separate_files=True)
def readHK(path):
    # N1 hk files are the 178-byte "hk_packet" defined in yingtian_packet.xml
    # (same 0x1a2b3c4d head magic as grid1x_hk_packet but different layout)
    return Dict(
        parse_grid_data_new(
            path,
            xml_file=str(Path(__file__).parent / "yingtian_packet.xml"),
            data_tag="hk_packet",
            endian="MSB",
        )[0]
    )


def getHK(sciFile, mode="ft"):
    """pair `{temp}-{bias...}-{idx}.event.dat` with `{temp}-{bias...}-ecu_*.hk`;
    when several candidates exist (e.g. CLYC m10C-265 has 3), pick the one
    whose utc range overlaps the sci utc range the most"""
    sciFile = Path(sciFile)
    parts = sciFile.stem.split("-")
    for n in (3, 2):
        prefix = "-".join(parts[:n])
        cands = [
            f
            for f in sciFile.parent.glob(f"{prefix}-ecu_*.hk")
            if "CI" not in f.stem
        ]
        if cands:
            break
    if not cands:
        raise FileNotFoundError(f"HK file for {sciFile} does not exist.")
    if len(cands) == 1:
        return cands[0]
    sci = readSci(str(sciFile), mode=mode)
    lo, hi = float(np.min(sci.utc)), float(np.max(sci.utc))
    best, best_overlap = None, -1.0
    for f in cands:
        tel = readHK(str(f))
        u = np.asarray(tel.utc_time, dtype=float)
        overlap = max(0.0, min(hi, u.max()) - max(lo, u.min()))
        if overlap > best_overlap:
            best, best_overlap = f, overlap
    return best


def _split_segments(utc):
    """split event utc into segments at the midpoints of safe gaps (no events
    within +-2 s of the cut), as in the gridN_cali L0 notebooks; segment count
    is adaptive (scan files can carry extra re-measured bias points)"""
    utc_sorted = np.sort(np.asarray(utc, dtype=float))
    dt = np.diff(utc_sorted)
    safe = np.where(dt > 4.0)[0]
    cuts = (utc_sorted[safe] + utc_sorted[safe + 1]) / 2.0
    return np.digitize(np.asarray(utc, dtype=float), cuts)


def single_read_n1(path: str, mode="ft", seg_bias=None, **kwargs):
    """read one N1 file; for bias-scan files pass seg_bias (bias code, e.g.
    285) to select the matching bias segment (mean bias within 0.1 V)."""
    sciFile = Path(path)
    hkFile = getHK(sciFile, mode)
    ow = kwargs.get("overwrite_cache", False)
    sci = readSci(str(sciFile), mode=mode, overwrite_cache=ow)
    tel = readHK(str(hkFile), overwrite_cache=ow)

    utc = np.asarray(sci.utc, dtype=float)
    # cut hk to the sci time coverage
    hk_utc = np.asarray(tel.utc_time, dtype=float)
    hk_in = (hk_utc >= utc.min()) & (hk_utc <= utc.max())

    if seg_bias is not None:
        seg = _split_segments(utc)
        n_seg = seg.max() + 1
        # mean hk bias (4-channel mean of vMon) per segment; select the
        # segment closest to the target bias
        b4 = np.stack([np.asarray(tel[f"sipm_voltage{i}"], dtype=float) / 1000
                       for i in range(4)]).mean(axis=0)
        codes = []
        for i in range(n_seg):
            lo_t, hi_t = utc[seg == i].min(), utc[seg == i].max()
            in_seg = (hk_utc >= lo_t) & (hk_utc <= hi_t)
            codes.append(float(np.mean(b4[in_seg])) if in_seg.any() else np.nan)
        codes = np.array(codes)
        target = int(np.nanargmin(np.abs(codes - seg_bias / 10.0)))
        if not np.isfinite(codes[target]) or abs(codes[target] - seg_bias / 10.0) > 0.1:
            raise ValueError(
                f"{path}: no segment at bias {seg_bias / 10.0} V "
                f"(segment means: {np.round(codes, 2).tolist()})"
            )
        sel = seg == target
        hk_in &= (hk_utc >= utc[sel].min()) & (hk_utc <= utc[sel].max())
    else:
        sel = np.ones(len(utc), dtype=bool)

    for k, v in list(sci.items()):
        v = np.asarray(v)
        if v.shape[0] == len(utc):
            sci[k] = v[sel]

    tel.tempSipm = [np.asarray(tel[f"sipm_temp{i}"]) / 100 - 273.15 for i in range(4)]
    tel.iMon = [np.asarray(tel[f"sipm_current{i}"]) for i in range(4)]
    tel.vMon = [np.asarray(tel[f"sipm_voltage{i}"]) / 1000 for i in range(4)]
    tel.bias = [tel.vMon[i] - 499 * tel.iMon[i] * 1e-6 for i in range(4)]
    for k in tel.keys():
        v = tel[k]
        if isinstance(v, list) and len(v) == 4:
            tel[k] = [np.asarray(v[i])[hk_in] for i in range(4)]

    if len(sci.data_max) == len(sci.data_base):
        sci.amp = sci.data_max - sci.data_base / 4.0
    else:
        assert False, f"{path} data_max.len != data_base.len"
    sci["timestampEvt"] = sci.timestamp

    n = sci.data_max.shape[0]
    for k, v in sci.items():
        if v.shape[0] == n and k != "channel_n":
            sci[k] = [v[sci.channel_n == i] for i in range(4)]
    return sci, tel
