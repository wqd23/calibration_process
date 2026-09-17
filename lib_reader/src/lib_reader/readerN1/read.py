# -*- coding:utf-8 -*-
"""Reader for the GRID-N1 (detect-Ecu.603) payload.

Two acquisition modes on the two crystal datasets:

- ``ft`` (GAGG): 584-byte ``grid1x_ft_packet`` -- 38 events/packet with a
  14-byte per-event stride (an extra ``data_ccm`` field).
- ``wf`` (CLYC): 1080-byte ``grid1x_wf_packet`` (512-sample waveform, defined in
  ``grid_packet_512wave.xml``), parsed in chunks with CRC-bad events dropped.

Both carry ``data_max``/``data_base``/``data_sum``, so
``amp = data_max - data_base/4`` (the integral charge ``quantity="q"`` is
``data_sum - data_base/4 * sample_length``).  HK is the 178-byte ``hk_packet``
from ``yingtian_packet.xml``.  N1 science packets carry a real ``utc``, so sci
events are matched to HK by nearest utc coverage and bias-scan files are split
into per-bias segments at the utc gaps; pass ``seg_bias`` (bias code, e.g.
``285``) to select the segment whose mean bias matches within 0.1 V.

Caches are rooted at ``data/<version>/`` where ``<version>`` is taken from the
raw path, so each N1 sub-version keeps its own L1/L2.
"""
from pathlib import Path

import numpy as np
from addict import Dict

from ..packet_parser import parse_grid_data_new, parse_grid_data_iter
from ..l1_cache import get_l1_frames, get_l2_processed, apply_selection

_SCI_XML = str(Path(__file__).with_name("grid_packet.xml"))
_WF_XML = str(Path(__file__).with_name("grid_packet_512wave.xml"))
_HK_XML = str(Path(__file__).with_name("yingtian_packet.xml"))


def _cache_ver(path) -> str:
    """Version component of ``.../data/<ver>/<root>/...`` (falls back to 'N1').

    ``<ver>`` may itself contain slashes (e.g. ``GRIDN1/GAGG``), so take every
    segment between ``data/`` and the data root (``raw_data``/``ec_src``/
    ``ec_xray``).
    """
    parts = Path(path).parts
    if "data" in parts:
        i = parts.index("data")
        roots = {"raw_data", "ec_src", "ec_xray"}
        segs = parts[i + 1:]
        for j, seg in enumerate(segs):
            if seg in roots:
                return "/".join(segs[:j]) or "N1"
    return "N1"


def _readSci_impl(path, mode="ft"):
    if mode == "ft":
        return Dict(parse_grid_data_new(
            path, xml_file=_SCI_XML, data_tag="grid1x_ft_packet", endian="MSB",
            multi_evt=38, multi_step=14)[0])
    if mode == "wf256":
        # the standard 568-byte 256-sample waveform packet (grid_packet.xml);
        # used by the neutron-beam runs
        return Dict(parse_grid_data_new(
            path, xml_file=_SCI_XML, data_tag="grid1x_wf_packet", endian="MSB")[0])
    if mode == "wf":
        parts = []
        for chunk in parse_grid_data_iter(
            path, xml_file=_WF_XML, data_tag="grid1x_wf_packet", endian="MSB",
            packet_len=1080, chunk_events=50000,
        ):
            d = chunk["data"]
            ok = np.asarray(d["crc_check"], dtype=bool)
            # TB/EC only use the onboard scalar sums; drop the 512-sample
            # waveform per chunk (keeping it across a multi-GB run would OOM)
            parts.append({k: np.asarray(v)[ok] for k, v in d.items()
                          if k != "waveform_data"
                          and np.asarray(v).shape[:1] == ok.shape[:1]})
        return Dict({k: np.concatenate([p[k] for p in parts]) for k in parts[0]})
    raise ValueError("Invalid mode. Use 'ft' or 'wf'.")


def _readHK_impl(path):
    return Dict(parse_grid_data_new(
        path, xml_file=_HK_XML, data_tag="hk_packet", endian="MSB")[0])


def _readSci(path, mode, overwrite, ver):
    return get_l1_frames(
        ver, "n1", path, {"sci": lambda: _readSci_impl(path, mode)},
        {"mode": mode}, overwrite=overwrite)["sci"]


def _readHK(path, overwrite, ver):
    return get_l1_frames(
        ver, "n1", path, {"hk": lambda: _readHK_impl(path)}, {},
        overwrite=overwrite)["hk"]


def getHK(sciFile, mode="ft", ver=None):
    """Pair ``{temp}-{bias...}-{idx}.event.dat`` with ``{temp}-{bias...}-ecu_*.hk``.

    When several candidates exist, pick the one whose utc range overlaps the
    science utc range the most.
    """
    sciFile = Path(sciFile)
    ver = ver or _cache_ver(str(sciFile))
    parts = sciFile.stem.split("-")
    cands = []
    for n in (3, 2):
        prefix = "-".join(parts[:n])
        cands = [f for f in sciFile.parent.glob(f"{prefix}-ecu_*.hk")
                 if "CI" not in f.stem]
        if cands:
            break
    if not cands:
        raise FileNotFoundError(f"HK file for {sciFile} does not exist.")
    if len(cands) == 1:
        return cands[0]
    sci = _readSci(str(sciFile), mode, False, ver)
    lo, hi = float(np.min(sci.utc)), float(np.max(sci.utc))
    best, best_overlap = None, -1.0
    for f in cands:
        tel = _readHK(str(f), False, ver)
        u = np.asarray(tel.utc_time, dtype=float)
        overlap = max(0.0, min(hi, u.max()) - max(lo, u.min()))
        if overlap > best_overlap:
            best, best_overlap = f, overlap
    return best


def _split_segments(utc):
    """Split event utc into segments at the midpoints of safe gaps (no events
    within +-2 s of the cut), as in the gridN_cali L0 notebooks."""
    utc_sorted = np.sort(np.asarray(utc, dtype=float))
    dt = np.diff(utc_sorted)
    safe = np.where(dt > 4.0)[0]
    cuts = (utc_sorted[safe] + utc_sorted[safe + 1]) / 2.0
    return np.digitize(np.asarray(utc, dtype=float), cuts)


def _single_readN1_impl(path, mode, hk_name, seg_bias, quantity, overwrite,
                        select=None, ver=None):
    ver = ver or _cache_ver(path)
    sci = _readSci(path, mode, overwrite, ver)
    if select is not None:
        sci = apply_selection(ver, "n1", path, {"mode": mode},
                              select[0], select[1], sci)
    tel = _readHK(str(hk_name), overwrite, ver)

    utc = np.asarray(sci.utc, dtype=float)
    hk_utc = np.asarray(tel.utc_time, dtype=float)
    hk_in = (hk_utc >= utc.min()) & (hk_utc <= utc.max())

    if seg_bias is not None:
        seg = _split_segments(utc)
        n_seg = seg.max() + 1
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
    # the HK monitor reads the regulated SiPM-side voltage (equals the file-name
    # setpoint), so no 499 ohm series-resistor drop is subtracted; pending
    # hardware confirmation, see docs/intermediate_data.md section 6.1
    tel.vMon = [np.asarray(tel[f"sipm_voltage{i}"]) / 1000 for i in range(4)]
    tel.bias = [tel.vMon[i] for i in range(4)]
    for k in list(tel.keys()):
        v = tel[k]
        if isinstance(v, list) and len(v) == 4:
            tel[k] = [np.asarray(v[i])[hk_in] for i in range(4)]

    if len(sci.data_max) == len(sci.data_base):
        if quantity == "q":
            if mode == "wf" and "sample_length" in sci:
                sample_length = np.asarray(sci.sample_length, dtype=float)
            else:
                sample_length = 256.0
            sci.amp = sci.data_sum - (sci.data_base / 4.0) * sample_length
        else:
            sci.amp = sci.data_max - sci.data_base / 4.0
    else:
        raise ValueError(f"{path}: data_max.len != data_base.len")
    sci["timestampEvt"] = sci.timestamp

    sci.pop("waveform_data", None)
    n = sci.data_max.shape[0]
    for k, v in list(sci.items()):
        v = np.asarray(v)
        if v.ndim == 1 and v.shape[0] == n and k != "channel_n":
            sci[k] = [v[sci.channel_n == i] for i in range(4)]
    return sci, tel


def single_readN1(path: str, config=None, mode="ft", seg_bias=None,
                  quantity="amp", hk_path=None, cache_ver=None, **kwargs):
    overwrite = kwargs.get("overwrite_cache", False)
    select = kwargs.get("select")
    ver = cache_ver or _cache_ver(path)
    params = {
        "mode": mode,
        "seg_bias": seg_bias,
        "quantity": quantity,
        "hk_path": None if hk_path is None else str(hk_path),
    }
    if select is not None:
        params["select"] = select[0]

    def process():
        hk_name = Path(hk_path) if hk_path else getHK(path, mode, ver)
        return _single_readN1_impl(path, mode, hk_name, seg_bias, quantity,
                                   overwrite, select=select, ver=ver)

    return get_l2_processed(ver, "n1", path, params, process,
                            overwrite=overwrite)
