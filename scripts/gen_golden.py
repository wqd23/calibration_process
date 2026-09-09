# -*- coding:utf-8 -*-
"""Generate the git golden files under ``tests/golden/`` (offline, one-time).

The committed golden is **self-contained**: it needs neither ``raw_data`` nor
``.oracle`` at test time.  A test reconstructs ``TBPoint`` / ``ECPoint`` lists
from the committed points JSON, runs the global-fit stage, and asserts the
produced coefficients equal the committed golden coefficients.

Points are faithful to the real pipeline: TB points come from the frozen
``.oracle`` single-fit pickles (they are dill-loadable), EC points come from a
freshly run single-fit store (the legacy EC pickles reference the removed
``operation`` module and are not loadable in the new package).  Golden
*coefficients* are always the authoritative legacy ``.oracle`` outputs.

Usage (data-dependent, run once; afterwards the committed golden is enough):

    # 1. populate a one-time single-fit store per EC representative version:
    #    python -m calibration_process.cli fit <ver> ec-src -o <store>/<ver>
    #    python -m calibration_process.cli fit <ver> ec-xray -o <store>/<ver>
    #    (TB points are read straight from .oracle, no store needed)
    # 2. python scripts/gen_golden.py [--store-dir /path/to/store]

Default store dir is /tmp/opencode/store; append the version.  See the
``STORE``/``TB_VERSIONS``/``EC_VERSIONS`` tables below for the version matrix.
"""

import argparse
import json
import os
import pickle
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src" / "calibration_process"))

from types import SimpleNamespace  # noqa: E402

from calibration_process import manifest as man  # noqa: E402
from calibration_process.pipeline import load_rt  # noqa: E402
from calibration_process.workflows import common as stages  # noqa: E402

GOLDEN = ROOT / "tests" / "golden"
ORACLE = ROOT / ".oracle"
CONFIG = ROOT / "src" / "calibration_process" / "configs"

# version -> which global-fit golden to generate.  Each representative exercises
# a distinct orchestration behaviour whose result is matched back to the frozen
# legacy oracle:
#   09  tb curvefit + ec polyfit(4ch)
#   11B tb lmfit + ec channel_count 3 (polyfit)
#   12B tb bias_min_filter
#   03B ec lmfit resolution
#   (04 exprfit-resolution and 10B channel_count-3 are NOT golden-tested here:
#    the legacy run produced no EC global output for them, so there is no
#    authoritative reference — pinning them would be self-referential.  They are
#    covered elsewhere: resolution methods + channel-count cap are unit-tested
#    in tests/test_stages.py.)
# One-time single-fit store root per version (EC points source; legacy EC
# pickles reference the removed ``operation`` module and cannot be loaded).
# ``--store-dir`` overrides the base directory (default /tmp/opencode/store);
# each version's store is ``<store_dir>/<ver>``.
STORE = {
    "09": Path("/tmp/opencode/all09_1788894811"),
    "03B": Path("/tmp/opencode/store/03B"),
    "11B": Path("/tmp/opencode/store/11B"),
}
STORE_BASE = Path("/tmp/opencode/store")
TB_VERSIONS = ["09", "11B", "12B"]
EC_VERSIONS = ["09", "03B", "11B"]


def _load_fp(path):
    d = pickle.load(open(path, "rb"))
    return SimpleNamespace(fit_result=d["fit_result"], tel=d["tel"])


def _latest(path, suffix):
    """Pick the highest timestamp (14-digit prefix) file matching suffix."""
    cands = [f for f in os.listdir(path) if f.endswith(suffix)]
    if not cands:
        return None
    return max(cands)


def _dump_points(points_by_channel, out_path):
    flat = [p for ch in points_by_channel for p in ch]
    data = [
        {
            "measurement_id": p.measurement_id,
            "channel": p.channel,
            "temperature": p.temperature,
            "temperature_err": p.temperature_err,
            "bias": p.bias,
            "bias_err": p.bias_err,
            "peak_center": p.peak_center,
            "peak_center_err": p.peak_center_err,
            "enabled": p.enabled,
        }
        for p in flat
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(data, f, indent=1)


def _dump_ec_points(src_pts, x_pts, out_path):
    flat = []
    for ch in src_pts:
        for p in ch:
            flat.append({"source_kind": "src", **ec_fields(p)})
    for ch in x_pts:
        for p in ch:
            flat.append({"source_kind": "xray", **ec_fields(p)})
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(flat, f, indent=1)


def ec_fields(p):
    return {
        "measurement_id": p.measurement_id,
        "channel": p.channel,
        "energy": p.energy,
        "peak_center": p.peak_center,
        "peak_center_err": p.peak_center_err,
        "resolution": p.resolution,
        "resolution_err": p.resolution_err,
        "enabled": p.enabled,
    }


def gen_tb(ver):
    rt = load_rt(ver)
    manifest = man.load_manifest(CONFIG / ver / "tb_manifest.yaml")
    items = [
        (m, _load_fp(ORACLE / ver / "TB_fit_result" / f"{Path(m.id).stem}.pickle"))
        for m in man.filtered_measurements(manifest)
    ]
    per_channel = stages.build_tb_points(rt, items)
    gdir = GOLDEN / ver
    _dump_points(per_channel, gdir / "points_tb.json")
    src = _latest(ORACLE / ver / "tb_logs", "_temp_bias_fit.json")
    if src:
        shutil.copy(ORACLE / ver / "tb_logs" / src, gdir / "tb_coeff.json")
    print(f"[{ver}] tb -> {gdir} ({len(per_channel[0])} points ch0)")


def gen_ec(ver, store_dir=None):
    rt = load_rt(ver)
    store = STORE[ver] if ver in STORE else (store_dir or STORE_BASE) / ver
    if not store.exists():
        raise FileNotFoundError(
            f"store for {ver} not found at {store}; populate it with\n"
            f"  python -m calibration_process.cli fit {ver} ec-src -o {store}\n"
            f"  python -m calibration_process.cli fit {ver} ec-xray -o {store}"
        )
    src_items, x_items = [], []
    for branch in ("ec_source", "ec_xray"):
        manifest = man.load_manifest(CONFIG / ver / f"{branch}_manifest.yaml")
        for m in man.filtered_measurements(manifest):
            fp = _load_fp(store / "single_process/EC_fit_result" / f"{Path(m.id).stem}.pickle")
            (src_items if branch == "ec_source" else x_items).append((m, fp))
    src_pts = stages.build_ec_points(rt, src_items, "src")
    x_pts = stages.build_ec_points(rt, x_items, "xray")
    gdir = GOLDEN / ver
    _dump_ec_points(src_pts, x_pts, gdir / "points_ec.json")
    elog = ORACLE / ver / "ec_logs"
    for ch in range(rt.payload.ec.channel_count):
        c = _latest(elog, f"_ec_coef_sci_ch{ch}.json")
        if c:
            shutil.copy(elog / c, gdir / f"ec_coeff_ch{ch}.json")
        d = _latest(elog, f"_ec_data_ch{ch}.npy")
        if d:
            shutil.copy(elog / d, gdir / f"ec_data_ch{ch}.npy")
    print(f"[{ver}] ec -> {gdir} (src={sum(len(c) for c in src_pts)} xray={sum(len(c) for c in x_pts)})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--store-dir", type=Path, default=STORE_BASE)
    args = parser.parse_args()
    for ver in TB_VERSIONS:
        gen_tb(ver)
    for ver in EC_VERSIONS:
        gen_ec(ver, store_dir=args.store_dir)


if __name__ == "__main__":
    main()
