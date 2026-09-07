# -*- coding:utf-8 -*-
"""High-level explicit pipeline for the new workflow.

These functions wire the shared stages together for a version/branch.  No
directory scanning happens here: measurement sets come from the manifest, and
every scientific computation delegates to the shared stages (which reuse the
protected kernel).
"""

from __future__ import annotations

from pathlib import Path

from . import manifest as man
from .runtime import load_runtime
from .workflows import common as stages

import os

CONFIG_ROOT = Path("src/calibration_process/configs")
DATA_ROOT = Path("data")


def _config_root(ver: str) -> Path:
    return CONFIG_ROOT / ver


def _manifest_path(ver: str, branch: str) -> Path:
    return _config_root(ver) / f"{branch}_manifest.yaml"


def output_root(ver: str, override: Path | str | None = None) -> Path:
    if override is not None:
        return Path(override)
    env = os.environ.get("CALIB_OUTPUT_ROOT")
    return Path(env) / ver if env else DATA_ROOT / ver


def branch_for_manifest(branch: str) -> str:
    return {
        "tb": "tb",
        "ec": "ec",
        "ec-src": "ec_source", "ecsrc": "ec_source", "ec_src": "ec_source",
        "ec_source": "ec_source",
        "ec-xray": "ec_xray", "ecxray": "ec_xray", "ec_xray": "ec_xray",
    }.get(branch, branch)


def load_rt(ver: str, out: Path | None = None):
    return load_runtime(ver, _config_root(ver), DATA_ROOT / ver, output_root(ver, out))


def discover(ver: str, branch: str) -> int:
    rt = load_rt(ver)
    raw_branch = branch_for_manifest(branch)
    draft = man.discover(ver, raw_branch, rt, DATA_ROOT / ver)
    man.dump_manifest(draft, _manifest_path(ver, raw_branch))
    return len(draft.measurements)


def list_measurements(ver: str, branch: str) -> list:
    return [m.id for m in man.load_manifest(_manifest_path(ver, branch_for_manifest(branch))).measurements]


def single_fit(ver: str, branch: str, measurement_id: str, nocache: bool = False,
               out: Path | None = None):
    rt = load_rt(ver, out)
    raw_branch = branch_for_manifest(branch)
    manifest = man.load_manifest(_manifest_path(ver, raw_branch))
    for m in man.filtered_measurements(manifest):
        if m.id == measurement_id:
            return stages.run_single_fit(rt, raw_branch, m, rt.output_root, nocache=nocache)
    raise KeyError(f"measurement {measurement_id!r} not in {ver}/{raw_branch} manifest")


def fit_branch(ver: str, branch: str, nocache: bool = False, out: Path | None = None) -> int:
    rt = load_rt(ver, out)
    raw_branch = branch_for_manifest(branch)
    manifest = man.load_manifest(_manifest_path(ver, raw_branch))
    n = 0
    for m in man.filtered_measurements(manifest):
        stages.run_single_fit(rt, raw_branch, m, rt.output_root, nocache=nocache)
        n += 1
    return n


def global_tb(ver: str, out: Path | None = None) -> list:
    rt = load_rt(ver, out)
    manifest = man.load_manifest(_manifest_path(ver, "tb"))
    items = []
    for m in man.filtered_measurements(manifest):
        fp = stages.load_single_fp_from_store(rt, "tb", m, rt.output_root)
        items.append((m, fp))
    per_channel = stages.build_tb_points(rt, items)
    return stages.global_tb(rt, per_channel, rt.output_root / "tb_logs")


def global_ec(ver: str, out: Path | None = None) -> list:
    rt = load_rt(ver, out)
    src_items, x_items = [], []
    for branch in ("ec_source", "ec_xray"):
        manifest = man.load_manifest(_manifest_path(ver, branch))
        for m in man.filtered_measurements(manifest):
            fp = stages.load_single_fp_from_store(rt, branch, m, rt.output_root)
            items = src_items if branch == "ec_source" else x_items
            items.append((m, fp))
    src_pts = stages.build_ec_points(rt, src_items, "src")
    x_pts = stages.build_ec_points(rt, x_items, "xray")
    return stages.global_ec(rt, src_pts, x_pts, rt.output_root / "ec_logs")


def all_version(ver: str, nocache: bool = False, out: Path | None = None) -> None:
    fit_branch(ver, "tb", nocache, out)
    global_tb(ver, out)
    fit_branch(ver, "ec-src", nocache, out)
    fit_branch(ver, "ec-xray", nocache, out)
    global_ec(ver, out)
