# -*- coding:utf-8 -*-
"""Typed intermediate products for the explicit workflow.

These are the stable boundaries between workflow stages.  They are derived
from, and never replace, the scientific kernel outputs (pickles / numpy
arrays).  Only the fields actually consumed by downstream stages are carried.
"""

from __future__ import annotations

from dataclasses import dataclass


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
