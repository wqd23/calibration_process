# -*- coding:utf-8 -*-
"""Strict Pydantic schemas for the new YAML configuration layer.

These schemas only validate the configuration layer.  They never re-model
scientific arrays or reader output.  Unknown fields are an error, so a typo
in a config file fails before the workflow starts.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, ConfigDict, field_validator


class _StrictBase(BaseModel):
    model_config = ConfigDict(extra="forbid")


class FitRangeSet(_StrictBase):
    """One fit-range file: measurement_id -> [ [lo,hi] x4 channels ].

    Each measurement must define exactly 4 channels of 2-element ranges, the
    same shape the legacy ``fit_range.json`` used (list of 4 ``[lo, hi]``).
    """

    measurements: Dict[str, List[List[float]]]

    @field_validator("measurements")
    @classmethod
    def _check_shape(cls, v: Dict[str, List[List[float]]]) -> Dict[str, List[List[float]]]:
        for mid, ranges in v.items():
            if len(ranges) != 4:
                raise ValueError(f"measurement {mid!r} must define exactly 4 channels, got {len(ranges)}")
            for ch, r in enumerate(ranges):
                if len(r) != 2:
                    raise ValueError(f"measurement {mid!r} channel {ch} range must be [lo, hi], got {r}")
                lo, hi = r
                if hi <= lo:
                    raise ValueError(f"measurement {mid!r} channel {ch} range [lo, hi] not ordered: {r}")
        return v


class BranchAnalysis(_StrictBase):
    """Per-branch human-editable analysis settings (background/model)."""

    default_bkg: Optional[str] = "lin"
    # measurement_id -> bkg_form override (None means no background)
    background_overrides: Dict[str, Optional[str]] = {}
    # measurement_id -> peak_form (only used by the 10B+ peak_form path)
    peak_overrides: Dict[str, str] = {}
    qa: Dict[str, float] = {}


class AnalysisSchema(_StrictBase):
    """The analysis.yaml file: background/model defaults and overrides."""

    tb: BranchAnalysis
    ec_source: BranchAnalysis
    ec_xray: BranchAnalysis


class TBParams(_StrictBase):
    reader: str
    bin_width: int
    adc_max: float
    science_dir: str
    fit_range_file: str
    # version-specific temp-bias global fit parameters
    tb_fit_p0: Optional[List[float]] = None
    tb_fit_maxfev: int = 10000
    bias_min_filter: Optional[float] = None


class ECParams(_StrictBase):
    reader: str
    bin_width: int
    adc_max: float
    x_path: str
    src_path: str
    energy_split_low: float = 49.0
    energy_split_high: float = 55.0
    ref_temp: float = 25.0
    ref_bias: float = 28.5
    tb_ref_path: str
    fit_range_file: str
    # measurement_id -> energy (keV)
    energy_map: Dict[str, float] = {}
    resolution_method: str = "polyfit"  # polyfit | exprfit | lmfit
    channel_count: int = 4
    # optional per-version extras
    src_bkg_map: Optional[Dict[str, List[str]]] = None
    xray_drop_energies: Optional[List[str]] = None
    xray_require_hk: bool = False
    xray_require_fit_range: bool = False
    xray_drop_old: bool = False
    xray_bkg_rotation: str = "circle"  # circle (ch+1) | fixed [1,2,0,0]


class PayloadSchema(_StrictBase):
    """The payload.yaml file: version-level scientific / reader parameters."""

    version: str
    tb: TBParams
    ec: ECParams


class ManifestEntry(_StrictBase):
    """One human-confirmed measurement in a manifest."""

    id: str
    branch: str
    science_files: List[str] = []
    hk_files: List[str] = []
    aux_files: List[str] = []
    metadata: Dict[str, Any] = {}
    use: bool = True
    channels: Optional[Dict[str, Dict[str, Any]]] = None


class ManifestSchema(_StrictBase):
    """A branch manifest: the human-confirmed list of measurements."""

    version: str
    branch: str
    measurements: List[ManifestEntry]

    @field_validator("measurements")
    @classmethod
    def _no_duplicate_ids(cls, v: List[ManifestEntry]) -> List[ManifestEntry]:
        seen = set()
        for m in v:
            if m.id in seen:
                raise ValueError(f"duplicated measurement id {m.id!r}")
            seen.add(m.id)
        return v
