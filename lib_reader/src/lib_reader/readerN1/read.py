# -*- coding:utf-8 -*-
"""Reader for the GRID-N1 (detect-Ecu.603) payload.

Science ``.event.dat`` files carry grid1x waveform (256 samples) or feature
packets; both include ``data_sum``/``data_ccm``.  ``amp = data_max -
data_base/4`` (note the factor 4, unlike the older A-family readers).  HK
files use the ``yingtian`` ``hk_packet`` with ``utc_time`` and per-channel
``sipm_voltage``/``sipm_current``/``sipm_temp``.

Caches are rooted at ``data/<version>/`` where ``<version>`` is taken from the
raw path, so each N1 sub-version (N1-Neutron, N1-Gamma-*) keeps its own L1/L2.
"""
from pathlib import Path

import numpy as np
from addict import Dict

from ..packet_parser import parse_grid_data_new
from ..l1_cache import get_l1_frames, get_l2_processed, apply_selection

_SCI_XML = str(Path(__file__).with_name("grid_packet.xml"))
_HK_XML = str(Path(__file__).with_name("yingtian_packet.xml"))


def _cache_ver(path) -> str:
    """Version component of ``.../data/<ver>/...`` (falls back to 'N1')."""
    parts = Path(path).parts
    if "data" in parts:
        i = parts.index("data")
        if i + 1 < len(parts):
            return parts[i + 1]
    return "N1"


def _readSci_impl(path, mode="wf"):
    if mode == "wf":
        return Dict(parse_grid_data_new(
            path, xml_file=_SCI_XML, data_tag="grid1x_wf_packet", endian="MSB")[0])
    if mode == "ft":
        return Dict(parse_grid_data_new(
            path, xml_file=_SCI_XML, data_tag="grid1x_ft_packet", endian="MSB",
            multi_evt=38, multi_step=14)[0])
    raise ValueError("Invalid mode. Use 'wf' or 'ft'.")


def _readHK_impl(path):
    return Dict(parse_grid_data_new(
        path, xml_file=_HK_XML, data_tag="hk_packet", endian="MSB")[0])


def _readSci(path, mode, overwrite):
    ver = _cache_ver(path)
    return get_l1_frames(
        ver, "n1", path, {"sci": lambda: _readSci_impl(path, mode)},
        {"mode": mode}, overwrite=overwrite)["sci"]


def _readHK(path, overwrite):
    ver = _cache_ver(path)
    return get_l1_frames(
        ver, "n1", path, {"hk": lambda: _readHK_impl(path)}, {},
        overwrite=overwrite)["hk"]


def getHK(sciFile):
    """Best-effort HK pairing; manifests should pass ``hk_path`` explicitly."""
    sciFile = Path(sciFile)
    same = sciFile.with_suffix(".hk")
    if same.exists():
        return same
    cands = sorted(sciFile.parent.glob("*ecu_*.hk"))
    if len(cands) == 1:
        return cands[0]
    raise FileNotFoundError(f"cannot pair HK for {sciFile}; set hk_path")


def _single_readN1_impl(path, mode, hk_name, sci_half, overwrite, select=None):
    sci = _readSci(path, mode, overwrite)
    if select is not None:
        sci = apply_selection(_cache_ver(path), "n1", path, {"mode": mode},
                              select[0], select[1], sci)
    if sci_half is not None:
        n_all = sci.data_max.shape[0]
        half = n_all // 2
        sl = slice(0, half) if sci_half == "first" else slice(half, n_all)
        for k, v in sci.items():
            if np.asarray(v).shape[0] == n_all:
                sci[k] = np.asarray(v)[sl]

    hk = _readHK(str(hk_name), overwrite)

    sci.amp = sci.data_max - sci.data_base / 4.0
    sci["timestampEvt"] = sci.timestamp

    tel = hk
    tel.tempSipm = [tel[f"sipm_temp{i}"] / 100 - 273.15 for i in range(4)]
    tel.iMon = [tel[f"sipm_current{i}"] for i in range(4)]
    tel.vMon = [tel[f"sipm_voltage{i}"] / 1000 for i in range(4)]
    tel.bias = list(tel.vMon)

    # HK can include ramp-up/warm-up records; keep the dominant 0.1 V x 1 C epoch
    b4 = np.stack([np.asarray(b) for b in tel.bias])
    t4 = np.stack([np.asarray(t) for t in tel.tempSipm])
    bm, tm = b4.mean(axis=0), t4.mean(axis=0)
    bb, tb = np.round(bm, 1), np.round(tm)
    keys, counts = np.unique(np.stack([bb, tb], axis=1), axis=0, return_counts=True)
    b0, t0 = keys[counts.argmax()]
    stable = (np.abs(bb - b0) < 0.05) & (np.abs(tb - t0) < 0.6)
    stable &= np.logical_and.reduce(np.abs(b4 - bm) < 0.3, axis=0)
    if np.any(stable):
        for k in list(tel.keys()):
            v = tel[k]
            if isinstance(v, list) and len(v) == 4:
                tel[k] = [np.asarray(v[i])[stable] for i in range(4)]
            elif isinstance(v, np.ndarray) and v.shape[0] == stable.shape[0]:
                tel[k] = v[stable]

    sci.pop("waveform_data", None)
    n = sci.data_max.shape[0]
    for k, v in sci.items():
        if hasattr(v, "shape") and np.ndim(v) == 1 and v.shape[0] == n and k != "channel_n":
            sci[k] = [v[sci.channel_n == i] for i in range(4)]
    return sci, tel


def single_readN1(path: str, config=None, mode="wf", hk_path=None, sci_half=None, **kwargs):
    overwrite = kwargs.get("overwrite_cache", False)
    select = kwargs.get("select")
    params = {
        "mode": mode,
        "hk_path": None if hk_path is None else str(hk_path),
        "sci_half": sci_half,
    }
    if select is not None:
        params["select"] = select[0]

    def process():
        hk_name = Path(hk_path) if hk_path else getHK(path)
        return _single_readN1_impl(path, mode, hk_name, sci_half, overwrite, select=select)

    return get_l2_processed(_cache_ver(path), "n1", path, params, process, overwrite=overwrite)
