# -*- coding:utf-8 -*-
"""Version -> workflow module registry.

Version selection happens exactly once here.  After a workflow module is
imported, no ``if version == ...`` branch should be needed downstream.
"""

from __future__ import annotations

import importlib

_MODULE = {
    "03B": "calibration_process.workflows.versions.v03B",
    "04": "calibration_process.workflows.versions.v04",
    "05B": "calibration_process.workflows.versions.v05B",
    "07": "calibration_process.workflows.versions.v07",
    "10B": "calibration_process.workflows.versions.v10B",
    "11B": "calibration_process.workflows.versions.v11B",
    "09": "calibration_process.workflows.versions.v09",
    "12B": "calibration_process.workflows.versions.v12B",
}


def get_workflow(version: str):
    if version not in _MODULE:
        raise KeyError(f"no explicit workflow for version {version!r}; "
                       f"available: {sorted(_MODULE)}")
    return importlib.import_module(_MODULE[version])
