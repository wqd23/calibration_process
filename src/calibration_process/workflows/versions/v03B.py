# -*- coding:utf-8 -*-
"""Explicit workflow for version 03B.

- TB: ``rundata`` files in the temp-bias dir, excluding baseline / CI / 50C;
  reader ``03b`` (binary + scienceConfig), bin_width 6, adc_max 16384.
- EC source: 3 hardcoded sources (Na22/Am241/Cs137) + hardcoded backgrounds,
  reader ``03b-src``, bin_width 10, lmfit resolution.
- EC X-ray: per-channel ``rundata`` files grouped by tube energy (name token 1),
  excluding CI; ``circle`` background rotation, lmfit resolution, split 49/51.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

SRC_LIST = [
    "src_Na22_20m_10cm_rundata2021-05-05-15-29-56.dat",
    "src_Am241_5m_10cm_rundata2021-05-05-15-15-48.dat",
    "src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat",
]
SRC_BKG = [
    "src_bkg_5m_10cm_rundata2021-05-05-15-54-31.dat",
    "src_bkg_5m_10cm_rundata2021-05-05-14-55-58.dat",
    "src_bkg_5m_10cm_rundata2021-05-05-12-33-57.dat",
]


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version == "03B", "v03B workflow used for non-03B version"
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
    files = [f for f in os.listdir(full)
             if "rundata" in f and "baseline" not in f and "CI" not in f and "50C" not in f]
    return [_rec(sci_dir, "tb", f) for f in files]


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    src_dir = rt.payload.ec.src_path
    return [_rec(src_dir, "ec_source", f, aux=[src_dir + "/" + b],
                 energy=rt.energies.get(f))
            for f, b in zip(SRC_LIST, SRC_BKG)]


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    x_ch = [f for f in os.listdir(full) if "_rundata" in f and "CI" not in f]
    x_list = sorted({f.split("_")[1] for f in x_ch})
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
