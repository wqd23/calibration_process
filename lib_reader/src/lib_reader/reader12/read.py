# -*- coding:utf-8 -*-
"""
Reader for GRID 12B payload data.

Differences from reader11:
- 12B .dat files contain grid1x_ft_packet (feature) packets, so the default
  read mode is 'ft' instead of 'wf'
- 12B HK files use the 187-byte grid1x_hk_packet (11B uses the 139-byte
  hk_grid1x_packet); the matching grid_packet.xml is shipped in this package
- file naming is `TB_m30C_270.dat` / `TB_m30C_270.hk` (setting in the name,
  actual values read from HK) or `{idx}_{kv}_ch{n}.dat` for X-ray files,
  so HK pairing is by identical stem with a (kv, channel) fallback
- ground tests have no GPS: sci utc is all zero, so no utc-based HK time cut.
  Instead, HK files include the bias ramp-up phase (and sometimes the warm-up
  after the point, e.g. ecu_1_088.hk runs from -16 C back to +29 C at the same
  bias), so HK records are cut to the dominant stable (bias, temperature)
  epoch: the most populous 0.1 V x 1 C bin of the per-record channel means.
  Pass `hk_bias` to instead select the epoch at a given target bias (needed
  when one HK file covers several bias points, e.g. ecu_1_093.hk)
- pass `hk_path` to use an explicit HK file instead of getHK() pairing, and
  `sci_half` ('first'/'second') to keep only half of the events in file order
  (needed when one observe file covers two bias points, e.g. 033_observe.dat)
"""
from pathlib import Path

import numpy as np
from addict import Dict

from ..packet_parser import parse_grid_data_new
from ..l1_cache import get_l1_frames, get_l2_processed, apply_selection

_XML = str(Path(__file__).with_name("grid_packet.xml"))


def _cache_ver(path) -> str:
    """Version component of ``.../data/<ver>/<root>/...`` (falls back to '12B').

    The same reader serves sibling payloads that share the 12B format (13B), so
    each keeps its own L1/L2 cache under ``data/<ver>/``.
    """
    parts = Path(path).parts
    if "data" in parts:
        i = parts.index("data")
        roots = {"raw_data", "ec_src", "ec_xray"}
        segs = parts[i + 1:]
        for j, seg in enumerate(segs):
            if seg in roots:
                return "/".join(segs[:j]) or "12B"
    return "12B"


def _readSci_impl(path, mode="ft"):
    if mode == "wf":
        return Dict(parse_grid_data_new(path, xml_file=_XML, data_tag="grid1x_wf_packet", endian="MSB")[0])
    elif mode == "ft":
        return Dict(
            parse_grid_data_new(
                path, xml_file=_XML, data_tag="grid1x_ft_packet", endian="MSB", multi_evt=41, multi_step=12
            )[0]
        )
    else:
        raise ValueError("Invalid mode. Use 'wf' or 'ft'.")


def _readHK_impl(path):
    return Dict(parse_grid_data_new(path, xml_file=_XML, data_tag="grid1x_hk_packet", endian="MSB")[0])


def readSci(path, mode="ft", overwrite_cache=False):
    return get_l1_frames(
        _cache_ver(path), "12b", path,
        {"sci": lambda: _readSci_impl(path, mode=mode)}, {"mode": mode},
        overwrite=overwrite_cache,
    )["sci"]


def readHK(path, overwrite_cache=False):
    return get_l1_frames(
        _cache_ver(path), "12b", path,
        {"hk": lambda: _readHK_impl(path)}, {},
        overwrite=overwrite_cache,
    )["hk"]


def getHK(sciFile):
    sciFile = Path(sciFile)
    # TB / source data: hk shares the stem of the dat file
    hkFile = sciFile.with_suffix(".hk")
    if hkFile.exists():
        return hkFile
    # X-ray data: {idx}_{kv}_ch{n}.dat pairs with {idx'}_{kv}_ch{n}.hk
    parts = sciFile.stem.split("_")
    if len(parts) >= 3:
        kv, ch = parts[1], parts[2]
        cands = [f for f in sciFile.parent.glob(f"*_{kv}_{ch}.hk") if "old" not in f.stem]
        if cands:
            # retakes exist for a few points; use the latest one
            cands.sort(key=lambda f: int(f.stem.split("_")[0]))
            return cands[-1]
        # no hk at this (kv, ch): caller skips this energy point
    raise FileNotFoundError(f"HK file for {sciFile} does not exist.")


def _single_read12_impl(path, mode, hk_name, hk_bias, sci_half, overwrite, select=None,
                        ver="12B"):
    sciExtracted = readSci(path, mode=mode, overwrite_cache=overwrite)
    if select is not None:
        sciExtracted = apply_selection(ver, "12b", path, {"mode": mode}, select[0], select[1], sciExtracted)
    telExtracted = readHK(str(hk_name), overwrite_cache=overwrite)

    # keep only one half of the events in file order (file order == time
    # order; the 32-bit timestamp wraps on long runs, so index is the robust
    # axis). Used when one observe file covers two bias points.
    if sci_half is not None:
        n_all = sciExtracted.data_max.shape[0]
        half = n_all // 2
        sl = slice(0, half) if sci_half == "first" else slice(half, n_all)
        for k, v in sciExtracted.items():
            if np.asarray(v).shape[0] == n_all:
                sciExtracted[k] = np.asarray(v)[sl]

    # amp
    if len(sciExtracted.data_max) == len(sciExtracted.data_base):
        amp = sciExtracted.data_max - sciExtracted.data_base / 4.0
    else:
        assert False, f"{path} data_max.len != data_base.len"
    sciExtracted.amp = amp

    telExtracted.tempSipm = [telExtracted[f"sipm_temp{i}"] / 100 - 273.15 for i in range(4)]
    # current, unit uA
    telExtracted.iMon = [telExtracted[f"sipm_current{i}"] for i in range(4)]
    # bias monitor, unit V.  The HK monitor of this payload reads the
    # regulated SiPM-side voltage (it equals the file-name setpoint), so the
    # legacy 499 ohm series-resistor drop is not subtracted (double-counting).
    # Pending hardware confirmation, see docs/intermediate_data.md section 6.1.
    telExtracted.vMon = [telExtracted[f"sipm_voltage{i}"] / 1000 for i in range(4)]
    telExtracted.bias = [telExtracted.vMon[i] for i in range(4)]

    sciExtracted["timestampEvt"] = sciExtracted.timestamp

    # no utc cut: ground data have utc == 0. Instead cut HK to a stable
    # (bias, temp) epoch: the epoch at hk_bias if given, else the dominant
    # 0.1 V x 1 C bin (HK files can include ramp-up, warm-up and retake
    # epochs at the same bias but other temperatures)
    b4 = np.stack([np.asarray(b) for b in telExtracted.bias])      # (4, n)
    t4 = np.stack([np.asarray(t) for t in telExtracted.tempSipm])  # (4, n)
    bm, tm = b4.mean(axis=0), t4.mean(axis=0)
    if hk_bias is not None:
        stable = np.abs(bm - hk_bias) < 0.3
    else:
        bb, tb = np.round(bm, 1), np.round(tm)
        keys, counts = np.unique(np.stack([bb, tb], axis=1), axis=0, return_counts=True)
        b0, t0 = keys[counts.argmax()]
        stable = (np.abs(bb - b0) < 0.05) & (np.abs(tb - t0) < 0.6)
    stable &= np.logical_and.reduce(np.abs(b4 - bm) < 0.3, axis=0)
    if not np.any(stable):
        raise ValueError(f"{path}: no stable-bias HK records found")
    for k in telExtracted.keys():
        v = telExtracted[k]
        if isinstance(v, list) and len(v) == 4:
            telExtracted[k] = [np.asarray(v[i])[stable] for i in range(4)]
        elif isinstance(v, np.ndarray) and v.shape[0] == stable.shape[0]:
            telExtracted[k] = v[stable]

    n = sciExtracted.data_max.shape[0]
    for k, v in sciExtracted.items():
        if v.shape[0] == n and k != "channel_n":
            sciExtracted[k] = [v[sciExtracted.channel_n == i] for i in range(4)]
    del n

    sciExtracted.pop("waveform_data", None)
    return sciExtracted, telExtracted


def single_read12(path: str, config=None, mode="ft", hk_path=None, hk_bias=None,
                  sci_half=None, cache_ver=None, **kwargs):
    overwrite = kwargs.get("overwrite_cache", False)
    select = kwargs.get("select")
    ver = cache_ver or _cache_ver(path)
    params = {
        "mode": mode,
        "hk_path": None if hk_path is None else str(hk_path),
        "hk_bias": hk_bias,
        "sci_half": sci_half,
    }
    if select is not None:
        params["select"] = select[0]

    def process():
        hk_name = Path(hk_path) if hk_path else getHK(path)
        return _single_read12_impl(path, mode, hk_name, hk_bias, sci_half, overwrite,
                                   select=select, ver=ver)

    return get_l2_processed(ver, "12b", path, params, process, overwrite=overwrite)
