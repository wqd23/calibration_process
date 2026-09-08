# -*- coding:utf-8 -*-
"""Explicit workflow for version 10B.

- TB: ``observe`` files in ``raw_data/tb_data``, excluding ``50C_265``; reader
  ``10b``, bin_width 6, adc_max 16384.
- EC source: 4 hardcoded sources + backgrounds, reader ``10b``, bin_width 4.
- EC X-ray: per-channel ``observe_*_ch{n}.dat`` grouped by tube energy (name
  token 2); ``fixed`` [ch1,ch2,ch0,ch0] background rotation.
- EC is physically 3 channels (channel_count=3) padded to 4 downstream.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

SRC_LIST = [
    "073_observe_Cs137.dat",
    "085_observe_Na22.dat",
    "089_observe_Am241.dat",
    "077_observe_Co60.dat",
]
SRC_BKG = [
    "075_observe_Cs137_bkg.dat",
    "086_observe_Na22_bkg.dat",
    "090_observe_Am241_bkg.dat",
    "078_observe_Co60_bkg.dat",
]


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version == "10B", "v10B workflow used for non-10B version"
    if branch == "tb":
        return _enumerate_tb(rt, data_dir)
    if branch == "ec_source":
        return _enumerate_ec_source(rt, data_dir)
    if branch == "ec_xray":
        return _enumerate_ec_xray(rt, data_dir)
    raise ValueError(f"unknown branch {branch!r}")


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    sci_dir = rt.payload.tb.science_dir
    full = data_dir / sci_dir
    files = [f for f in os.listdir(full) if "observe" in f and "50C_265" not in f]
    return [_rec(sci_dir, "tb", f) for f in files]


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    src_dir = rt.payload.ec.src_path
    return [_rec(src_dir, "ec_source", f, aux=[src_dir + "/" + b],
                 energy=rt.energies.get(f))
            for f, b in zip(SRC_LIST, SRC_BKG)]


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    x_ch = [f for f in os.listdir(full) if "_ch" in f and "hk" not in f]
    x_list = sorted({f.split("_")[2] for f in x_ch})
    out = []
    for energy in x_list:
        ch_files = []
        for i in range(4):
            matched = [f for f in x_ch if f"{energy}_ch{i}" in f]
            ch_files.append(f"{x_dir}/{matched[0]}")
        out.append({"id": energy, "branch": "ec_xray", "science_files": ch_files,
                    "hk_files": [], "aux_files": [],
                    "metadata": {"energy": rt.energies.get(energy)}, "use": True})
    return out


def _rec(dir_, branch, name, aux=None, energy=None):
    return {"id": name, "branch": branch, "science_files": [f"{dir_}/{name}"],
            "hk_files": [], "aux_files": aux or [], "metadata": {"energy": energy},
            "use": True}
