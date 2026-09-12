# -*- coding:utf-8 -*-
"""Explicit workflow for the GRID-N1 sub-versions.

GRID-N1 is split by acquisition dataset into separate sub-version namespaces
that share this module:

- ``N1-Gamma-Am241``: the GAGG runs (Am241 59.5 keV, ft packets, ch1/2).
- ``N1-Gamma-Na22``: the CLYC runs (Na22 511 keV, wf 512-sample packets,
  ch0-3).
- ``N1-Gamma-EC`` / ``N1-Neutron``: placeholders for later phases.

TB measurement selection reproduces ``gridN_cali``'s curation (also used by the
legacy ``TB_operation_N1``): per dataset, temperature-bias runs are either
single-bias files (``{temp}-{bias}-{idx}.event.dat``) or bias-scan files
(``{temp}-265-290-{idx}.event.dat``) split into one segment per bias code via
the reader's ``seg_bias``.  CLYC scans carry extra re-measured bias segments.
Files with "To" (temperature transition), "CI" (current scan) or "test" are
excluded.  Measurement ids are ``{dataset}_{temp}_{bias}`` (matching the
committed ``fit_range.json``).
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np

from ...runtime import RuntimeConfig

VERSIONS = ("N1-Gamma-Am241", "N1-Gamma-Na22", "N1-Gamma-EC", "N1-Neutron")

# which dataset each gamma sub-version reads (TB)
DATASET = {"N1-Gamma-Am241": "GAGG", "N1-Gamma-Na22": "CLYC"}

# standard bias-scan codes and the CLYC extra re-measured segments, taken from
# the gridN_cali read_raw_CLYC.ipynb vol_sets
STANDARD_BIAS = (265, 270, 275, 280, 283, 285, 287, 290)
CLYC_EXTRA_BIAS = {
    "m20C": (277, 279, 281),
    "m10C": (279, 281),
    "0C": (281,),
    "10C": (289,),
    "20C": (289, 291),
    "30C": (289, 291),
}

# neutron/gamma PSD cut on ccm and amplitude (tune against real data)
NEUTRON_SELKEY = "n1-neutron@rev1"
NEUTRON_CCM_MIN = 0.0
NEUTRON_AMP_MIN = 0.0

_SINGLE = re.compile(r"^(m?\d+C)-(\d{3})-\d+\.event\.dat$")
_SCAN = re.compile(r"^(m?\d+C)-265-290-\d+\.event\.dat$")
# neutron-beam runs: {temperature}-{bias V}-{idx}.event.dat (a single bias)
_NEUTRON = re.compile(r"^(?P<temp>.+?)-(?P<bias>\d+\.\d+)-(?P<idx>\d+)\.event\.dat$")


def enumerate_measurements(version: str, branch: str, rt: RuntimeConfig,
                           data_dir: Path) -> list:
    assert version in VERSIONS, f"vN1 workflow used for non-N1 version {version!r}"
    if version == "N1-Neutron" and branch == "tb":
        return _enumerate_neutron_tb(rt, data_dir)
    if branch == "tb" and version in DATASET:
        return _enumerate_tb(rt, data_dir, DATASET[version])
    return []


def _enumerate_neutron_tb(rt: RuntimeConfig, data_dir: Path) -> list:
    """The neutron-beam temperature scans (fixed bias) used for a TB snapshot."""
    sci_dir = rt.payload.tb.science_dir
    full = Path(data_dir) / sci_dir
    out = []
    for f in sorted(os.listdir(full)):
        if not f.endswith(".event.dat") or "To" in f or "CI" in f or "test" in f:
            continue
        m = _NEUTRON.match(f)
        if not m:
            continue
        out.append({
            "id": f[: -len(".event.dat")],
            "branch": "tb",
            "science_files": [f"{sci_dir}/{f}"],
            "hk_files": [],
            "aux_files": [],
            "metadata": {"temp": m.group("temp"), "bias": float(m.group("bias"))},
            "use": True,
        })
    return out


def _enumerate_tb(rt: RuntimeConfig, data_dir: Path, dataset: str) -> list:
    sci_dir = rt.payload.tb.science_dir
    full = Path(data_dir) / sci_dir
    # keyed by id with last-wins, matching the legacy point_map (a directory can
    # hold both a single-bias file and a scan that cover the same point)
    point_map: dict = {}
    for f in sorted(os.listdir(full)):
        if not f.endswith(".event.dat") or "To" in f or "CI" in f or "test" in f:
            continue
        m = _SINGLE.match(f)
        if m:
            temp, bias = m.group(1), int(m.group(2))
            rec = _record(dataset, temp, bias, None, sci_dir, f)
            point_map[rec["id"]] = rec
            continue
        m = _SCAN.match(f)
        if m:
            temp = m.group(1)
            biases = list(STANDARD_BIAS)
            if dataset == "CLYC":
                biases += list(CLYC_EXTRA_BIAS.get(temp, ()))
            for bias in biases:
                rec = _record(dataset, temp, bias, bias, sci_dir, f)
                point_map[rec["id"]] = rec
    return list(point_map.values())


def _record(dataset, temp, bias, seg_bias, sci_dir, fname):
    return {
        "id": f"{dataset}_{temp}_{bias}",
        "branch": "tb",
        "science_files": [f"{sci_dir}/{fname}"],
        "hk_files": [],
        "aux_files": [],
        "metadata": {"dataset": dataset, "temp": temp, "bias": bias,
                     "seg_bias": seg_bias},
        "use": True,
    }


def selection(version: str, branch: str):
    """Version/branch event-selection hook (only the future N1-Neutron uses it)."""
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
