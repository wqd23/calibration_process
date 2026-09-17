# -*- coding:utf-8 -*-
"""Explicit workflow for the 15B payload (``cali_data/14B15B/02``).

TB files are ``{temp}_{bias}_{idx}.event.dat`` (mostly) or
``{temp}-{bias}-{idx}.event.dat`` (the m20C group), paired with
``{temp}_{bias}_ecu_*.hk`` / ``{temp}-{bias}-ecu_*.hk``.
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

VERSIONS = ("15B",)

_TB = re.compile(r"^(?P<temp>m?\d+C?)[-_](?P<bias>\d{3})[-_]\d+\.event\.dat$")


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> List[dict]:
    assert version in VERSIONS, f"v15B workflow used for non-15B version {version!r}"
    if branch == "tb":
        return _enumerate_tb(rt, data_dir)
    return []


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    sci_dir = rt.payload.tb.science_dir
    full = Path(data_dir) / sci_dir
    out = []
    for name in sorted(os.listdir(full)):
        m = _TB.match(name)
        if m is None:
            continue
        temp, bias = m.group("temp"), int(m.group("bias"))
        if not temp.endswith("C"):
            temp += "C"
        stem = Path(name).stem
        if "-" in stem:
            prefix = stem.rsplit("-", 1)[0]
        else:
            prefix = stem.rsplit("_", 1)[0]
        hks = sorted(full.glob(f"{prefix}-ecu_*.hk")) + sorted(full.glob(f"{prefix}_ecu_*.hk"))
        if not hks:
            continue
        out.append({
            "id": f"{temp}_{bias}",
            "branch": "tb",
            "science_files": [f"{sci_dir}/{name}"],
            "hk_files": [f"{sci_dir}/{hks[0].name}"],
            "aux_files": [],
            "metadata": {"hk_bias": bias / 10.0},
            "use": True,
        })
    return out
