# -*- coding:utf-8 -*-
"""
QA statistics: goodness-of-fit metrics distribution over processed results
----------
Scans single_process/TB_fit_result, single_process/EC_fit_result and
tb_logs of each payload version, summarizes redchi / success /
boundary_hit per channel, so that QA thresholds can be set from the
observed distribution of known-good results.

Usage:
    python3 -m calibration_process.qa_stats [ver ...]
"""

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from . import util_lib as util
from .__init__ import CFG_PATH

cfg = util.json_load(CFG_PATH)


def _single_fit_stats(root: Path):
    """Per-channel metrics from single spectrum fit pickles under root."""
    redchi = defaultdict(list)
    stats = {"ok": 0, "failed": 0, "no_range": 0, "boundary": 0}
    failed_files = []
    boundary_cases = []
    for pkl in sorted(root.glob("*.pickle")):
        try:
            data = util.pickle_load(pkl)
        except Exception:
            stats["failed"] += 1
            failed_files.append(f"{pkl.name} (unreadable)")
            continue
        fit_result = data.get("fit_result", [])
        for ich, res in enumerate(fit_result):
            if res is None:
                stats["no_range"] += 1
                continue
            if isinstance(res, dict) and res.get("success") is False:
                stats["failed"] += 1
                failed_files.append(f"{pkl.name} ch{ich}: {res.get('error', '')[:60]}")
                continue
            if not isinstance(res, dict) or "redchi" not in res:
                # legacy pickle without QA metrics
                continue
            stats["ok"] += 1
            redchi[ich].append(res["redchi"])
            if res["boundary_hit"]:
                stats["boundary"] += 1
                boundary_cases.append(f"{pkl.name} ch{ich}: {res['boundary_hit']}")
    return redchi, stats, failed_files, boundary_cases


def _percentile_str(values):
    if not values:
        return "n=%-4d (no data)" % 0
    q = np.percentile(values, [50, 90, 99, 100])
    return "n=%-4d p50=%.2f p90=%.2f p99=%.2f max=%.2f" % (
        len(values),
        q[0],
        q[1],
        q[2],
        q[3],
    )


def _tb_2d_stats(ver: str):
    """redchi of temp-bias 2D surface fits in tb_logs."""
    best = None
    for f in sorted(Path(f"data/{ver}/tb_logs").glob("*_temp_bias_fit.json")):
        best = f
    if best is None:
        return None
    redchis = [ch.get("redchi") for ch in util.json_load(str(best))]
    return best.name, redchis


def qa_stats(ver: str) -> None:
    print(f"\n==== {ver} ====")
    for kind, sub in (("tb", "TB_fit_result"), ("ec", "EC_fit_result")):
        redchi, stats, failed, boundary = _single_fit_stats(
            Path(f"data/{ver}/single_process/{sub}")
        )
        print(
            f"  [{kind}] ok={stats['ok']} failed={stats['failed']} "
            f"no_range={stats['no_range']} boundary_hit={stats['boundary']}"
        )
        for ich in sorted(redchi):
            print(f"    ch{ich} redchi: {_percentile_str(redchi[ich])}")
        for item in failed[:5]:
            print(f"    FAILED: {item}")
        if len(failed) > 5:
            print(f"    ... and {len(failed) - 5} more failures")
        for item in boundary[:3]:
            print(f"    BOUNDARY: {item}")
    tb2d = _tb_2d_stats(ver)
    if tb2d is not None:
        name, redchis = tb2d
        print(f"  [tb 2D] {name}: redchi per ch = {[round(r, 1) for r in redchis]}")


def main():
    vers = sys.argv[1:] or list(cfg.keys())
    for ver in vers:
        qa_stats(ver)


if __name__ == "__main__":
    main()
