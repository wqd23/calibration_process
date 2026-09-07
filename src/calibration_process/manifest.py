# -*- coding:utf-8 -*-
"""Manifest discovery and loading.

The scientific workflow only accepts a human-confirmed manifest.  Directory
scanning is only used for ``discover -> manifest draft``.  Runtime code never
re-scans directories to decide which measurements take part in analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import yaml

from .config_schema import ManifestEntry, ManifestSchema


# A raw measurement record as produced by a version's selection logic.
# shape: {id, branch, science_files, hk_files, aux_files, metadata, use, channels}
RawMeasurement = Dict


def version_data_dir(version: str) -> Path:
    return Path("data") / version


def _resolve_relative(data_dir: Path, rel: str) -> str:
    """Return the absolute path for a manifest-relative path.

    Paths stored in a manifest are relative to ``data/<version>``.  For 09
    the science files live under ``raw_data/...``, so the relative path e.g.
    ``raw_data/temp_bias/xxx.txt`` resolves to ``data/09/raw_data/...``.
    """
    p = data_dir / rel
    return str(p)


def discover(version: str, branch: str, payload, data_dir: Path) -> ManifestSchema:
    """Produce a manifest draft for a version/branch.

    The actual measurement enumeration is delegated to the version workflow
    module (``workflows.versions``), which encodes the historical selection
    rules.  This function only wraps the enumerated records into the schema
    and dedupes/validates ids.
    """
    from .workflows.registry import get_workflow

    wf = get_workflow(version)
    raw = wf.enumerate_measurements(version, branch, payload, data_dir)
    entries = []
    for r in raw:
        entry = ManifestEntry(
            id=r["id"],
            branch=r["branch"],
            science_files=list(r.get("science_files", [])),
            hk_files=list(r.get("hk_files", [])),
            aux_files=list(r.get("aux_files", [])),
            metadata=dict(r.get("metadata", {})),
            use=bool(r.get("use", True)),
            channels=r.get("channels"),
        )
        entries.append(entry)
    return ManifestSchema(version=version, branch=branch, measurements=entries)


def load_manifest(path: Path) -> ManifestSchema:
    return ManifestSchema.model_validate(yaml.safe_load(open(path)))


def dump_manifest(manifest: ManifestSchema, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(manifest.model_dump(), f, sort_keys=False)


def filtered_measurements(manifest: ManifestSchema) -> List[ManifestEntry]:
    """Return measurements with use=True."""
    return [m for m in manifest.measurements if m.use]
