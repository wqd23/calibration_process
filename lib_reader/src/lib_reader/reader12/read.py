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
from cachier import cachier

from .parse_grid_data import parse_grid_data_new


@cachier(cache_dir=Path(".cache") / "12B", separate_files=True)
def readSci(path, mode="ft"):
    if mode == "wf":
        return Dict(parse_grid_data_new(path, data_tag="grid1x_wf_packet", endian="MSB")[0])
    elif mode == "ft":
        return Dict(
            parse_grid_data_new(
                path, data_tag="grid1x_ft_packet", endian="MSB", multi_evt=41, multi_step=12
            )[0]
        )
    else:
        raise ValueError("Invalid mode. Use 'wf' or 'ft'.")


@cachier(cache_dir=Path(".cache") / "12B", separate_files=True)
def readHK(path):
    return Dict(parse_grid_data_new(path, data_tag="grid1x_hk_packet", endian="MSB")[0])


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


def single_read12(path: str, mode="ft", hk_path=None, hk_bias=None, sci_half=None, **kwargs):
    observe_name = path
    hk_name = Path(hk_path) if hk_path else getHK(observe_name)
    sciExtracted = readSci(observe_name, mode=mode, overwrite_cache=kwargs.get("overwrite_cache", False))
    telExtracted = readHK(str(hk_name), overwrite_cache=kwargs.get("overwrite_cache", False))

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
    # bias monitor, unit V
    telExtracted.vMon = [telExtracted[f"sipm_voltage{i}"] / 1000 for i in range(4)]
    telExtracted.bias = [
        telExtracted.vMon[i] - 499 * telExtracted.iMon[i] * 1e-6 for i in range(4)
    ]

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

    return sciExtracted, telExtracted
