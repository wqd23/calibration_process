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

from .config_schema import AnalysisSchema, FitRangeSet, PayloadSchema, ReaderSchema


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
    reader: Optional[ReaderSchema] = None
    # branch -> measurement_id -> list-of-4 [lo, hi]
    fit_ranges: Dict[str, Dict[str, List[List[float]]]] = field(default_factory=dict)
    # branch -> measurement_id -> bkg_form (resolved: override or default)
    bkg_forms: Dict[str, Dict[str, Optional[str]]] = field(default_factory=dict)
    # measurement_id -> energy (keV)
    energies: Dict[str, float] = field(default_factory=dict)
    # lazily built per-channel temp-bias correction functions (see ``corr``)
    _corr: Optional[List] = field(default=None, repr=False, compare=False)

    @property
    def corr(self) -> List:
        """Per-channel temp-bias gain correction, built on first access.

        Reading the ``tb_ref_path`` reference is deliberately deferred: a
        payload with no EC branch (or a TB-only sub-version whose reference is
        not produced yet) must still load its runtime.  Only the EC branches
        touch this, and a missing reference then fails with a clear message.
        """
        if self._corr is None:
            self._corr = _build_corr(self.payload, self.data_dir)
        return self._corr

    def fit_range(self, branch: str, measurement_id: str) -> List[List[float]]:
        return self.fit_ranges[branch][measurement_id]

    def bkg_form(self, branch: str, measurement_id: str) -> Optional[str]:
        return self.bkg_forms[branch].get(measurement_id, self.analysis_default(branch))

    def analysis_default(self, branch: str) -> Optional[str]:
        b = {
            "tb": self.analysis.tb,
            "ec_source": self.analysis.ec_source,
            "ec_xray": self.analysis.ec_xray,
        }.get(branch)
        return b.default_bkg if b is not None else None

    def reader_handler(self, ending: str) -> Optional[str]:
        """Registry name of the reader handler for ``ending`` (or None)."""
        if self.reader and ending in self.reader.readers:
            return self.reader.readers[ending].handler
        return None

    def reader_params(self, ending: str) -> Dict:
        """Version-specific reader parameters for ``ending`` (empty if none)."""
        if self.reader and ending in self.reader.readers:
            return dict(self.reader.readers[ending].params)
        return {}


def _resolve_bkg(branch: str, analysis: AnalysisSchema) -> Dict[str, Optional[str]]:
    b = {
        "tb": analysis.tb,
        "ec_source": analysis.ec_source,
        "ec_xray": analysis.ec_xray,
    }.get(branch)
    return dict(b.background_overrides) if b is not None else {}


def load_runtime(version: str, config_root: Path, data_dir: Path,
                 output_root: Optional[Path] = None) -> RuntimeConfig:
    payload = PayloadSchema.model_validate(
        _load_yaml(config_root / "payload.yaml")
    )
    analysis = AnalysisSchema.model_validate(
        _load_yaml(config_root / "analysis.yaml")
    )
    reader = None
    reader_path = config_root / "reader.yaml"
    if reader_path.exists():
        reader = ReaderSchema.model_validate(_load_yaml(reader_path))
    rt = RuntimeConfig(
        version=version,
        data_dir=data_dir,
        output_root=output_root or data_dir,
        config_root=config_root,
        payload=payload,
        analysis=analysis,
        reader=reader,
        energies=dict(payload.ec.energy_map),
    )
    # fit range files: tb, ec_source, ec_xray, neutron
    fit_files = {
        "tb": "fit_range_tb.yaml",
        "ec_source": "fit_range_ec_source.yaml",
        "ec_xray": "fit_range_ec_xray.yaml",
        "neutron": "fit_range_neutron.yaml",
    }
    for branch, fname in fit_files.items():
        p = config_root / fname
        if not p.exists():
            continue
        rt.fit_ranges[branch] = FitRangeSet.model_validate(
            _load_yaml(p)
        ).measurements
    for branch in ("tb", "ec_source", "ec_xray", "neutron"):
        rt.bkg_forms[branch] = _resolve_bkg(branch, analysis)
    # ``rt.corr`` is built lazily on first access (EC branches only), so a
    # TB-only / neutron payload does not need its TB reference to exist yet.
    return rt


def _build_corr(payload: PayloadSchema, data_dir: Path) -> List:
    from . import util_lib as util

    tb_ref = data_dir / payload.ec.tb_ref_path
    if not tb_ref.exists():
        raise util.FitError(
            f"EC correction reference not found: {tb_ref}; run `calib global "
            f"<ver> tb` (and merge per-source results where applicable) first"
        )
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
