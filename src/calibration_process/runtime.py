# -*- coding:utf-8 -*-
"""Resolved runtime configuration for one version.

Loads the YAML configs (payload / analysis / fit-ranges) once and exposes
resolved, version-independent structures.  Resolved config is never dumped to
disk.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import yaml

from .config_schema import AnalysisSchema, FitRangeSet, PayloadSchema


def _load_yaml(path: Path):
    with open(path) as f:
        return yaml.safe_load(f)


@dataclass
class RuntimeConfig:
    version: str
    data_dir: Path
    output_root: Path
    config_root: Path
    payload: PayloadSchema
    analysis: AnalysisSchema
    # branch -> measurement_id -> list-of-4 [lo, hi]
    fit_ranges: Dict[str, Dict[str, List[List[float]]]] = field(default_factory=dict)
    # branch -> measurement_id -> bkg_form (resolved: override or default)
    bkg_forms: Dict[str, Dict[str, Optional[str]]] = field(default_factory=dict)
    # measurement_id -> energy (keV)
    energies: Dict[str, float] = field(default_factory=dict)
    # per EC branch corr functions (ref temp-bias correction)
    corr: List = field(default_factory=list)

    def fit_range(self, branch: str, measurement_id: str) -> List[List[float]]:
        return self.fit_ranges[branch][measurement_id]

    def bkg_form(self, branch: str, measurement_id: str) -> Optional[str]:
        return self.bkg_forms[branch].get(measurement_id, self.analysis_default(branch))

    def analysis_default(self, branch: str) -> Optional[str]:
        b = {
            "tb": self.analysis.tb,
            "ec_source": self.analysis.ec_source,
            "ec_xray": self.analysis.ec_xray,
        }[branch]
        return b.default_bkg


def _resolve_bkg(branch: str, analysis: AnalysisSchema) -> Dict[str, Optional[str]]:
    b = {
        "tb": analysis.tb,
        "ec_source": analysis.ec_source,
        "ec_xray": analysis.ec_xray,
    }[branch]
    return dict(b.background_overrides)


def load_runtime(version: str, config_root: Path, data_dir: Path,
                 output_root: Optional[Path] = None) -> RuntimeConfig:
    payload = PayloadSchema.model_validate(
        _load_yaml(config_root / "payload.yaml")
    )
    analysis = AnalysisSchema.model_validate(
        _load_yaml(config_root / "analysis.yaml")
    )
    rt = RuntimeConfig(
        version=version,
        data_dir=data_dir,
        output_root=output_root or data_dir,
        config_root=config_root,
        payload=payload,
        analysis=analysis,
        energies=dict(payload.ec.energy_map),
    )
    # fit range files: order tb, ec_source, ec_xray
    fit_files = {
        "tb": "fit_range_tb.yaml",
        "ec_source": "fit_range_ec_source.yaml",
        "ec_xray": "fit_range_ec_xray.yaml",
    }
    for branch, fname in fit_files.items():
        p = config_root / fname
        if not p.exists():
            continue
        rt.fit_ranges[branch] = FitRangeSet.model_validate(
            _load_yaml(p)
        ).measurements
    for branch in ("tb", "ec_source", "ec_xray"):
        rt.bkg_forms[branch] = _resolve_bkg(branch, analysis)
    # corr from TB reference
    rt.corr = _build_corr(payload, data_dir)
    return rt


def _build_corr(payload: PayloadSchema, data_dir: Path) -> List:
    from . import util_lib as util

    tb_ref = data_dir / payload.ec.tb_ref_path
    tb_result = util.json_load(str(tb_ref))
    ref_temp = payload.ec.ref_temp
    ref_bias = payload.ec.ref_bias
    ref_func = [
        lambda t, b: util.tempbias2DFunction(
            t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
        )
        for c in tb_result
    ]
    return [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
