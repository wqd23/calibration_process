# -*- coding:utf-8 -*-
"""Explicit workflow for the 14B payload (``cali_data/14B15B/03``).

TB: one measurement per ``Temp_Vbias/{temp}-{bias}-{idx}.event.dat``, paired
with ``{temp}-{bias}-ecu_*.hk``.  The target bias comes from the file-name
code, so no extra metadata is needed.  Source/X-ray (EC) branches are not
enumerated yet.
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig

VERSIONS = ("14B",)

_TB = re.compile(r"^(m?\d+C)-(\d{3})-\d+\.event\.dat$")


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> List[dict]:
    assert version in VERSIONS, f"v14B workflow used for non-14B version {version!r}"
    if branch == "tb":
        return _enumerate_tb(rt, data_dir)
    return []


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    sci_dir = rt.payload.tb.science_dir
    full = Path(data_dir) / sci_dir
    out = []
    for name in sorted(os.listdir(full)):
        if not _TB.match(name):
            continue
        prefix = name.rsplit("-", 1)[0]
        hks = sorted(full.glob(f"{prefix}-ecu_*.hk"))
        out.append({
            "id": name[: -len(".event.dat")],
            "branch": "tb",
            "science_files": [f"{sci_dir}/{name}"],
            "hk_files": [f"{sci_dir}/{hks[-1].name}"] if hks else [],
            "aux_files": [],
            "metadata": {},
            "use": True,
        })
    return out
