# -*- coding:utf-8 -*-
"""Explicit workflow for version 09.

Reproduces the historical selection rules and pipeline for 09:

- TB: every ``.txt`` in ``raw_data/temp_bias`` except 6 explicitly excluded
  files; temperature-bias global fit with the default curvefit parameters.
- EC source: 3 hardcoded sources (Na22 / Co60 / Cs137) with hardcoded
  backgrounds; energy split at the 49/55 keV K edge.
- EC X-ray: one 4-channel file set per tube energy, excluding ``65keV_``;
  the 4 channels are read as separate files and reconstructed (fp03B path).
"""

from __future__ import annotations

import os
from pathlib import Path

from ...runtime import RuntimeConfig

# --- version-specific selection constants ---------------------------------- #
TB_EXCLUDE = [
    "0826_30C_265_5m_0x0090.txt",
    "0826_30C_265_5m_0x00A0.txt",
    "0826_30C_265_5m_0x00B0.txt",
    "0826_30C_265_5m_0x00BB.txt",
    "0826_20C_290_5m_0x01F0 (1).txt",
    "0825_0C_290_4m_0x00C5.txt",
]

SRC_LIST = [
    "0826_10C_285_Na22_20m_0x00CF.txt",
    "0827_10C_285_Co60_20m_0x010F.txt",
    "0827_10C_285_Cs137_20m_0x010F.txt",
]
SRC_BKG = [
    "0826_10C_285_Na22_bkg_20m.txt",
    "0827_10C_285_bkg_20m_0x00CF.txt",
    "0827_10C_285_bkg_20m_0x00CF.txt",
]


def _rel(prefix: str, name: str) -> str:
    return f"{prefix}/{name}"


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    """Return the manifest-raw measurement records for a branch of version 09."""
    assert version == "09", "v09 workflow used for non-09 version"
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
    for excl in TB_EXCLUDE:
        files.remove(excl)
    out = []
    for f in files:
        out.append({
            "id": f,
            "branch": "tb",
            "science_files": [_rel(sci_dir, f)],
            "hk_files": [],
            "aux_files": [],
            "metadata": {},
            "use": True,
        })
    return out


def _enumerate_ec_source(rt: RuntimeConfig, data_dir: Path) -> list:
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


def _enumerate_ec_xray(rt: RuntimeConfig, data_dir: Path) -> list:
    x_dir = rt.payload.ec.x_path
    full = data_dir / x_dir
    drop = set(rt.payload.ec.xray_drop_energies or [])
    x_ch = [f for f in os.listdir(full) if "_ch" in f and f.split("_")[0] not in drop]
    x_list = sorted({f.split("_")[0] for f in x_ch})
    out = []
    for energy_name in x_list:
        ch_files = []
        for i in range(4):
            matched = [f for f in x_ch if f"{energy_name}_ch{i}" in f]
            ch_files.append(_rel(x_dir, matched[0]))
        out.append({
            "id": energy_name,
            "branch": "ec_xray",
            "science_files": ch_files,
            "hk_files": [],
            "aux_files": [],
            "metadata": {"energy": rt.energies.get(energy_name)},
            "use": True,
        })
    return out
