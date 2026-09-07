# -*- coding:utf-8 -*-
"""Version -> workflow module registry.

Version selection happens exactly once here.  After a workflow module is
imported, no ``if version == ...`` branch should be needed downstream.
"""

from __future__ import annotations

import importlib

_MODULE = {
    "09": "calibration_process.workflows.versions.v09",
    # only 09 is fully migrated; add more versions as their Gate B completes
}


def get_workflow(version: str):
    if version not in _MODULE:
        raise KeyError(f"no explicit workflow for version {version!r}; "
                       f"available: {sorted(_MODULE)}")
    return importlib.import_module(_MODULE[version])
