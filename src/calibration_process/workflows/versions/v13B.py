# -*- coding:utf-8 -*-
"""Explicit workflow for the 13B payload (device 073).

TB files live in ``raw_data`` under two naming styles:

- ``{idx}_TB_{temp}_{bias}.dat`` with the HK at a different index
  (``{idx'}_TB_{temp}_{bias}.hk``), e.g. the 20C group;
- ``TB_{temp}_{bias}.dat`` / ``TB_{temp}_{bias}.hk``.

They are paired by ``(temp, bias)`` and the target bias from the code is passed
as ``hk_bias``.  A point whose HK is missing (e.g. 20C/270 only ships a log
``.txt``) is skipped.
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

VERSIONS = ("13B",)

_A = re.compile(r"^\d+_TB_(?P<temp>m?\d+C)_(?P<bias>\d{3})\.dat$")
_B = re.compile(r"^TB_(?P<temp>m?\d+C)_(?P<bias>\d{3})\.dat$")


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> List[dict]:
    assert version in VERSIONS, f"v13B workflow used for non-13B version {version!r}"
    if branch == "tb":
        return _enumerate_tb(rt, data_dir)
    return []


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    sci_dir = rt.payload.tb.science_dir
    full = Path(data_dir) / sci_dir
    out = []
    for name in sorted(os.listdir(full)):
        m = _A.match(name) or _B.match(name)
        if m is None:
            continue
        temp, bias = m.group("temp"), int(m.group("bias"))
        hks = sorted(full.glob(f"*TB_{temp}_{bias}.hk"))
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
