# -*- coding:utf-8 -*-
"""Explicit workflow for version 04.

- TB: every ``.txt`` in ``raw_data/20210501_tempbias_Am241_GRID04``.
- EC source: 4 hardcoded sources (Co60/Na22/Cs137/Am241) + hardcoded bkg.
- EC X-ray: per-channel ``.txt`` files grouped by tube energy (name token 3),
  excluding ``_18p0_``; 4-channel reconstruction (fp03B), ``circle`` bkg
  rotation, ``ExprFit`` resolution fit.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

SRC_LIST = [
    "210504162829_COM6_src_Co60_10m_10cm-Data.txt",
    "210504170604_COM6_src_Na22_30m_10cm-Data.txt",
    "210505120306_COM6_src_Cs137_15m_10cm-Data.txt",
    "210505151551_COM6_src_Am241_5m_10cm-Data.txt",
]
SRC_BKG = [
    "210504164227_COM6_src_bkg_5m_10cm-Data.txt",
    "210504175217_COM6_src_bkg_5m_10cm-Data.txt",
    "210505123401_COM6_src_bkg_5m_10cm-Data.txt",
    "210505145547_COM6_src_bkg_5m_10cm-Data.txt",
]


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version == "04", "v04 workflow used for non-04 version"
    if branch == "tb":
        return _enumerate_tb(rt, data_dir)
    if branch == "ec_source":
        return _enumerate_ec_source(rt, data_dir)
    if branch == "ec_xray":
        return _enumerate_ec_xray(rt, data_dir)
    raise ValueError(f"unknown branch {branch!r}")


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path) -> list:
    sci_dir = rt.payload.tb.science_dir
    full = data_dir / sci_dir
    files = [f for f in os.listdir(full) if os.path.splitext(f)[1] == ".txt"]
    out = []
    for f in files:
        out.append({"id": f, "branch": "tb",
                    "science_files": [f"{sci_dir}/{f}"],
                    "hk_files": [], "aux_files": [], "metadata": {}, "use": True})
    return out


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> list:
    src_dir = rt.payload.ec.src_path
    out = []
    for f, b in zip(SRC_LIST, SRC_BKG):
        out.append({"id": f, "branch": "ec_source",
                    "science_files": [f"{src_dir}/{f}"],
                    "hk_files": [], "aux_files": [f"{src_dir}/{b}"],
                    "metadata": {"energy": rt.energies.get(f)}, "use": True})
    return out


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    subs = set(rt.payload.ec.xray_drop_substrs or [])
    x_ch = [f for f in os.listdir(full) if "_ch" in f and not any(s in f for s in subs)]
    x_list = sorted({f.split("_")[3] for f in x_ch})
    out = []
    for energy in x_list:
        ch_files = []
        for i in range(rt.payload.ec.channel_count):
            matched = [f for f in x_ch if f"{energy}_ch{i}" in f]
            ch_files.append(f"{x_dir}/{matched[0]}")
        out.append({
            "id": energy, "branch": "ec_xray", "science_files": ch_files,
            "hk_files": [], "aux_files": [],
            "metadata": {"energy": rt.energies.get(energy)}, "use": True,
        })
    return out
