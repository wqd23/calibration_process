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


_BRANCHES = ("tb", "ec_source", "ec_xray")


def present_branches(ver: str) -> list:
    """The branches that have a manifest for ``ver`` (order: tb, ec_source, ec_xray).

    Payload versions need not provide every branch: a gamma-only or neutron
    sub-version simply omits the manifests it does not use.
    """
    return [b for b in _BRANCHES if _manifest_path(ver, b).exists()]


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


def build_tb_points(ver: str, out: Path | None = None):
    """L4: build per-channel TB points from the L3/L2 artefacts."""
    rt = load_rt(ver, out)
    manifest = man.load_manifest(_manifest_path(ver, "tb"))
    items = []
    for m in man.filtered_measurements(manifest):
        fp = stages.load_single_fp_from_store(rt, "tb", m, rt.output_root)
        items.append((m, fp))
    return rt, stages.build_tb_points(rt, items)


def build_ec_points(ver: str, out: Path | None = None):
    """L4: build per-channel EC points (source + xray) from the L3/L2 artefacts."""
    rt = load_rt(ver, out)
    src_items, x_items = [], []
    for branch in ("ec_source", "ec_xray"):
        mp = _manifest_path(ver, branch)
        if not mp.exists():
            continue
        manifest = man.load_manifest(mp)
        for m in man.filtered_measurements(manifest):
            fp = stages.load_single_fp_from_store(rt, branch, m, rt.output_root)
            items = src_items if branch == "ec_source" else x_items
            items.append((m, fp))
    src_pts = stages.build_ec_points(rt, src_items, "src")
    x_pts = stages.build_ec_points(rt, x_items, "xray")
    return rt, src_pts, x_pts


def global_tb(ver: str, out: Path | None = None) -> list:
    rt, per_channel = build_tb_points(ver, out)
    return stages.global_tb(rt, per_channel, rt.output_root / "tb_logs")


def global_ec(ver: str, out: Path | None = None) -> list:
    rt, src_pts, x_pts = build_ec_points(ver, out)
    return stages.global_ec(rt, src_pts, x_pts, rt.output_root / "ec_logs")


# Declarative step list: L1/L2 are produced together by the reader call (L1
# faithful frames + L2 processed output), L3 is the single fit, L4 builds the
# points and L5 runs the global fits.  ``--until`` stops after a given layer.
STEP_ORDER = ("L1", "L2", "L3", "L4", "L5")


def process_version(ver: str, nocache: bool = False, out: Path | None = None) -> None:
    """L1+L2 only: run the readers for every measurement without fitting."""
    branches = present_branches(ver)
    if not branches:
        raise FileNotFoundError(
            f"{ver}: no *_manifest.yaml found under {_config_root(ver)}"
        )
    rt = load_rt(ver, out)
    for branch in branches:
        manifest = man.load_manifest(_manifest_path(ver, branch))
        for m in man.filtered_measurements(manifest):
            fc = stages.single_run_spec(rt, branch, m)
            stages.build_fit_operation(rt, branch, fc, nocache=nocache)


def all_version(ver: str, nocache: bool = False, out: Path | None = None,
                until: str = "L5") -> None:
    until = until.upper()
    if until not in STEP_ORDER:
        raise ValueError(f"unknown layer {until!r}; choose from {STEP_ORDER}")
    branches = present_branches(ver)
    if not branches:
        raise FileNotFoundError(
            f"{ver}: no *_manifest.yaml found under {_config_root(ver)}"
        )
    idx = STEP_ORDER.index(until)
    if idx <= STEP_ORDER.index("L2"):
        process_version(ver, nocache, out)
        return
    for manifest_branch, cli_branch in (
        ("tb", "tb"), ("ec_source", "ec-src"), ("ec_xray", "ec-xray"),
    ):
        if manifest_branch in branches:
            fit_branch(ver, cli_branch, nocache, out)
    if idx == STEP_ORDER.index("L3"):
        return
    rt = load_rt(ver, out)
    tb_points = src_pts = x_pts = None
    if "tb" in branches:
        _rt_tb, tb_points = build_tb_points(ver, out)
    if "ec_source" in branches or "ec_xray" in branches:
        _rt_ec, src_pts, x_pts = build_ec_points(ver, out)
    if idx == STEP_ORDER.index("L4"):
        return
    if tb_points is not None:
        stages.global_tb(rt, tb_points, rt.output_root / "tb_logs")
    if src_pts is not None or x_pts is not None:
        stages.global_ec(rt, src_pts, x_pts, rt.output_root / "ec_logs")
