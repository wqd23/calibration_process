# -*- coding:utf-8 -*-
"""Explicit workflow for the GRID-N1 sub-versions.

GRID-N1 is split by source into four namespaces that share this module:

- ``N1-Gamma-Am241`` / ``N1-Gamma-Na22``: TB temperature-bias scans (the
  per-source TB results are merged by ``scripts/merge_tb.py`` for the EC ref).
- ``N1-Gamma-EC``: energy calibration against the merged TB reference.
- ``N1-Neutron``: a standalone peak-fit profile; neutron/gamma separation uses
  the packet ``data_ccm`` (PSD) selection hook below.

The raw measurement selection is manual (the gridN_cali reference curates its
L1 configs by hand), so ``enumerate_measurements`` reads a committed
``<branch>_file_map.json`` under the version's config directory.  Each entry::

    {"id": "...", "science_files": ["raw_data/..."],
     "hk_files": ["raw_data/..."], "aux_files": [], "metadata": {}}
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ...runtime import RuntimeConfig

VERSIONS = ("N1-Gamma-Am241", "N1-Gamma-Na22", "N1-Gamma-EC", "N1-Neutron")

# neutron/gamma PSD cut on ccm and amplitude (tune against real data)
NEUTRON_SELKEY = "n1-neutron@rev1"
NEUTRON_CCM_MIN = 0.0
NEUTRON_AMP_MIN = 0.0


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version in VERSIONS, f"vN1 workflow used for non-N1 version {version!r}"
    map_path = Path(rt.config_root) / f"{branch}_file_map.json"
    if not map_path.exists():
        return []
    entries = json.loads(map_path.read_text())
    out = []
    for e in entries:
        out.append({
            "id": e["id"],
            "branch": branch,
            "science_files": list(e.get("science_files", [])),
            "hk_files": list(e.get("hk_files", [])),
            "aux_files": list(e.get("aux_files", [])),
            "metadata": dict(e.get("metadata", {})),
            "use": bool(e.get("use", True)),
        })
    return out


def selection(version: str, branch: str):
    """Version/branch event-selection hook (only N1-Neutron selects)."""
    if version == "N1-Neutron" and branch == "neutron":
        return _neutron_selection()
    return None


def _neutron_selection():
    def fn(frames):
        amp = np.asarray(frames["data_max"]) - np.asarray(frames["data_base"]) / 4.0
        ccm = np.asarray(frames["data_ccm"]) / 65535.0
        mask = (amp >= NEUTRON_AMP_MIN) & (ccm >= NEUTRON_CCM_MIN)
        return mask, {"amp": amp, "ccm": ccm}

    return NEUTRON_SELKEY, fn
