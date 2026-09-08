# -*- coding:utf-8 -*-
"""Explicit workflow for version 12B.

Reproduces the historical selection rules for 12B:

- TB: measurement set comes from ``tb_file_map.json`` (backup dir), with two
  excluded points (-20,275)/(-20,285), a custom temp-bias initial guess / maxfev,
  and the global fit restricted to bias >= 27.25 V.  Two points share one
  observe file and are split by ``sci_half`` with a target ``hk_bias``.
- EC source: 4 sources, all sharing the environmental background 0611env.dat.
- EC X-ray: one 4-channel file set per tube kV (x_path), dropping ``old``
  retakes and any point without a complete 4-channel HK pairing or a complete
  4-channel fit range.  Background comes from the ``fixed`` [ch1,ch2,ch0,ch0]
  rotation.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

EXCLUDE = {(-20, 275), (-20, 285)}

SRC_LIST = [
    "0611_Na22_10min_240f0032.dat",
    "0611_Cs137_30min_240f0064.dat",
    "0611_Co60_25min_240f0032.dat",
    "0611_Am241_12min_240f0032.dat",
]
SRC_BKG = ["0611env.dat"] * 4


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> List[dict]:
    assert version == "12B", "v12B workflow used for non-12B version"
    if branch == "tb":
        return _enumerate_tb(rt, data_dir)
    if branch == "ec_source":
        return _enumerate_ec_source(rt, data_dir)
    if branch == "ec_xray":
        return _enumerate_ec_xray(rt, data_dir)
    raise ValueError(f"unknown branch {branch!r}")


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    file_map = json.load(open(data_dir / rt.payload.tb.tb_file_map))
    point_map = {}
    for p in file_map:
        if (p["temp_setpoint_C"], p["bias_code"]) in EXCLUDE:
            continue
        name = f"{p['temp_setpoint_C']}C_{p['bias_code']}"
        point_map[name] = p
    out = []
    for name in sorted(point_map):
        p = point_map[name]
        note = p.get("note", "")
        meta = {}
        if "前半段" in note:
            meta = {"sci_half": "first", "hk_bias": p["bias_setpoint_V"]}
        elif "后半段" in note:
            meta = {"sci_half": "second", "hk_bias": p["bias_setpoint_V"]}
        out.append({
            "id": name,
            "branch": "tb",
            "science_files": ["raw_data/" + p["observe_file"]],
            "hk_files": ["raw_data/" + p["hk_file"]],
            "aux_files": [],
            "metadata": meta,
            "use": True,
        })
    return out


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    src_dir = rt.payload.ec.src_path
    out = []
    for f, bkg in zip(SRC_LIST, SRC_BKG):
        out.append({
            "id": f,
            "branch": "ec_source",
            "science_files": [_rel(src_dir, f)],
            "hk_files": [],
            "aux_files": [_rel(src_dir, bkg)],
            "metadata": {"energy": rt.energies.get(f)},
            "use": True,
        })
    return out


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    drop_old = rt.payload.ec.xray_drop_old
    x_ch = [f for f in os.listdir(full)
            if f.endswith(".dat") and "_ch" in f and (not drop_old or "old" not in f)]
    x_list = sorted({f.split("_")[1] for f in x_ch}, key=int)
    out = []
    for energy in x_list:
        if rt.payload.ec.xray_require_hk and not _x_hk_complete(full, x_ch, energy):
            continue
        if rt.payload.ec.xray_require_fit_range and not _range_complete(rt, energy):
            continue
        ch_files = [_rel(x_dir, _get_x_file(x_ch, energy, i)) for i in range(4)]
        out.append({
            "id": energy,
            "branch": "ec_xray",
            "science_files": ch_files,
            "hk_files": [],
            "aux_files": [],
            "metadata": {"energy": rt.energies.get(energy)},
            "use": True,
        })
    return out


def _get_x_file(x_ch, energy: str, i: int) -> str:
    matched = [f for f in x_ch if f.split("_")[1] == energy and f"_ch{i}" in f]
    return matched[0]


def _x_hk_complete(x_dir: Path, x_ch, energy: str) -> bool:
    from lib_reader.reader12.read import getHK

    try:
        for i in range(4):
            getHK(str(x_dir / _get_x_file(x_ch, energy, i)))
    except (IndexError, FileNotFoundError):
        return False
    return True


def _range_complete(rt: RuntimeConfig, energy: str) -> bool:
    ranges = rt.fit_ranges.get("ec_xray", {}).get(energy)
    if ranges is None:
        return False
    return all(r is not None for r in ranges)


def _rel(prefix: str, name: str) -> str:
    return f"{prefix}/{name}"
