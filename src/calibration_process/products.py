# -*- coding:utf-8 -*-
"""Typed intermediate products for the explicit workflow.

These are the stable boundaries between workflow stages.  They are derived
from, and never replace, the scientific kernel outputs (pickles / numpy
arrays).  Only the fields actually consumed by downstream stages are carried.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Union


@dataclass(frozen=True)
class FileRunSpec:
    """Resolved read/spectrum/fit config for one single-file run.

    Mirrors the legacy ``file_config()`` tuple
    ``(read_config, bkg_read_config, spectrum_config, fit_config)``.
    """

    # file_lib.Read_config
    read_config: object
    # file_lib.Read_config
    bkg_read_config: object
    # file_lib.Spectrum_config
    spectrum_config: object
    # file_lib.Fit_config
    fit_config: object


@dataclass(frozen=True)
class SingleFitResult:
    """Per-channel result of a single peak fit (derived from the pickle)."""

    measurement_id: str
    channel: int
    peak_amplitude: float
    peak_amplitude_err: float
    peak_center: float
    peak_center_err: float
    peak_sigma: float
    peak_sigma_err: float
    resolution: float
    resolution_err: float
    rate: float
    rate_err: float
    redchi: float
    ndf: int
    success: bool
    qa_flag: str
    boundary_hit: List[str]


@dataclass(frozen=True)
class TBPoint:
    """One temperature-bias point for one channel, for TB global fit."""

    measurement_id: str
    channel: int
    temperature: float
    temperature_err: float
    bias: float
    bias_err: float
    peak_center: float
    peak_center_err: float
    enabled: bool = True


@dataclass(frozen=True)
class ECPoint:
    """One energy-calibration point for one channel, for EC global fit."""

    measurement_id: str
    channel: int
    source_kind: str  # "xray" | "src"
    energy: float
    peak_center: float
    peak_center_err: float
    resolution: float
    resolution_err: float
    enabled: bool = True


@dataclass(frozen=True)
class MeasurementBundle:
    """A human-confirmed measurement record (from a manifest)."""

    id: str
    branch: str  # "tb" | "ec_source" | "ec_xray"
    science_files: List[str]  # manifest-relative paths
    hk_files: List[str]
    aux_files: List[str]
    metadata: dict
    use: bool = True
    channels: dict = None  # channel -> {"use": bool}

    def channel_use(self, channel: int) -> bool:
        if not self.use:
            return False
        if self.channels and str(channel) in self.channels:
            return bool(self.channels[str(channel)].get("use", True))
        return True


OptionalFit = Union[SingleFitResult, None]
FitResultList = List[OptionalFit]
