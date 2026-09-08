# -*- coding:utf-8 -*-
"""Explicit workflow for version 07.

- TB: every ``.txt`` in ``raw_data/北师大正样温度偏压标定数据-20211124``.
- EC source: sources matched by ``src`` / ``_bk_`` / ``bkg`` name filters, with
  a per-source background discovered from ``bkg{src}``.
- EC X-ray: per-channel `.txt` grouped by tube energy (name token 2), dropping
  the ``40p0/15p0/12p0/99p9/90p1`` points; ``circle`` rotation, polyfit
  resolution.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

XRAY_DROP = ["40p0", "15p0", "12p0", "99p9", "90p1"]


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version == "07", "v07 workflow used for non-07 version"
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
    files = [f for f in os.listdir(full) if os.path.splitext(f)[1] == ".txt"]
    out = []
    for f in files:
        out.append({"id": f, "branch": "tb", "science_files": [f"{sci_dir}/{f}"],
                    "hk_files": [], "aux_files": [], "metadata": {}, "use": True})
    return out


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    src_dir = rt.payload.ec.src_path
    full = data_dir / src_dir
    src_list = [f for f in os.listdir(full)
                if "src" in f and "_bk_" not in f and "bkg" not in f]
    out = []
    for f in src_list:
        src = f.split("_")[3][:2]
        bkg = [g for g in os.listdir(full)
               if "_bk_" not in g and f"bkg{src}" in g][0]
        out.append({"id": f, "branch": "ec_source",
                    "science_files": [f"{src_dir}/{f}"],
                    "hk_files": [], "aux_files": [f"{src_dir}/{bkg}"],
                    "metadata": {"energy": rt.energies.get(f)}, "use": True})
    return out


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    x_ch = [f for f in os.listdir(full)
            if "_ch" in f and not any(s in f for s in XRAY_DROP)]
    x_list = sorted({f.split("_")[2] for f in x_ch})
    out = []
    for energy in x_list:
        ch_files = []
        for i in range(rt.payload.ec.channel_count):
            matched = [f for f in x_ch if energy in f and f"_ch{i}_" in f]
            ch_files.append(f"{x_dir}/{matched[0]}")
        out.append({
            "id": energy, "branch": "ec_xray", "science_files": ch_files,
            "hk_files": [], "aux_files": [],
            "metadata": {"energy": rt.energies.get(energy)}, "use": True,
        })
    return out
