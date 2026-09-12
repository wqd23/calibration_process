# -*- coding:utf-8 -*-
"""Populate the L1/L2 parquet caches for every measurement (read step only).

For each version/branch it resolves the run spec and builds the
``File_operation_05b`` (whose constructor performs the data readout via
``read_out``), which populates ``data/{ver}/l1/`` (faithful frames) and
``data/{ver}/l2/`` (processed ``(sci, tel)``) without running any fit.
Errors on individual measurements are reported and skipped.

Usage:
    python scripts/warm_l1_cache.py [ver ...]
"""
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src" / "calibration_process"))

from calibration_process import manifest as man  # noqa: E402
from calibration_process.pipeline import load_rt  # noqa: E402
from calibration_process.workflows import common as stages  # noqa: E402

CONFIG = ROOT / "src" / "calibration_process" / "configs"
BRANCHES = ["tb", "ec_source", "ec_xray", "neutron"]
ALL_VERSIONS = ["03B", "04", "05B", "07", "09", "10B", "11B", "12B",
                "GRIDN1/GAGG", "GRIDN1/CLYC", "GRIDN1/EC", "GRIDN1/Neutron"]


def warm_version(ver: str):
    rt = load_rt(ver)
    n_ok = n_skip = n_fail = 0
    for branch in BRANCHES:
        mpath = CONFIG / ver / f"{branch}_manifest.yaml"
        if not mpath.exists():
            continue
        manifest = man.load_manifest(mpath)
        for m in man.filtered_measurements(manifest):
            try:
                fc = stages.single_run_spec(rt, branch, m)
                stages.build_fit_operation(rt, branch, fc)
                n_ok += 1
            except FileNotFoundError:
                n_skip += 1
            except Exception as e:  # noqa: BLE001
                n_fail += 1
                print(f"  FAIL {ver}/{branch}/{m.id}: {type(e).__name__}: {e}")
                traceback.print_exc()
    print(f"[{ver}] read ok={n_ok} skip={n_skip} fail={n_fail}")


def main():
    versions = sys.argv[1:] or ALL_VERSIONS
    for ver in versions:
        warm_version(ver)


if __name__ == "__main__":
    main()
