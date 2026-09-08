# -*- coding:utf-8 -*-
"""Explicit workflow for version 05B.

- TB: ``rundata`` files in the temperature-bias dir, excluding 50C.
- EC source: ``rundata`` source files (no ``bkg`` in the name), each paired to
  a ``src_bkg`` file sharing the source token.
- EC X-ray: **one 4-channel file per tube setting** (not per-channel files);
  distinguished by ``_observe.dat``, excluding ``XM_22``.  Background is the
  same file read with a cyclically shifted per-channel time cut, and the
  reader is the ``xray`` reader driven by ``x_config``.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List

from ...runtime import RuntimeConfig


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version == "05B", "v05B workflow used for non-05B version"
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
    files = [f for f in os.listdir(full) if "rundata" in f and "50C" not in f]
    return [_rec(sci_dir, "tb", f) for f in files]


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    src_dir = rt.payload.ec.src_path
    full = data_dir / src_dir
    src_list = [f for f in os.listdir(full) if "rundata" in f and "bkg" not in f]
    out = []
    for f in src_list:
        token = f.split("_")[1]
        bkg = [g for g in os.listdir(full) if token in g and "bkg" in g and "rundata" in g][0]
        out.append({"id": f, "branch": "ec_source",
                    "science_files": [f"{src_dir}/{f}"], "hk_files": [],
                    "aux_files": [f"{src_dir}/{bkg}"],
                    "metadata": {"energy": rt.energies.get(f)}, "use": True})
    return out


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> List[dict]:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    files = [f for f in os.listdir(full) if "_observe.dat" in f and "XM_22" not in f]
    return [_rec(x_dir, "ec_xray", f, energy=rt.energies.get(f)) for f in files]


def _rec(dir_, branch, name, energy=None):
    return {"id": name, "branch": branch,
            "science_files": [f"{dir_}/{name}"], "hk_files": [], "aux_files": [],
            "metadata": {"energy": energy}, "use": True}
