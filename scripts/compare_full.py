# -*- coding:utf-8 -*-
"""Full legacy vs new 09 comparison (Gate B evidence).

Compares single-fit pickles, tb/ec global outputs (json/npy), and figure
file sets + dimensions between the frozen legacy oracle (.oracle/09) and the
new workflow output (data/09).
"""

import os
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "tests")
from compare import assert_pickle_equivalent, Mismatch  # noqa: E402

import sys as _sys
VER = _sys.argv[1] if len(_sys.argv) > 1 else "09"
ORACLE = Path(f".oracle/{VER}")
NEW = Path(f".new/{VER}")


def strip_ts(name: str) -> str:
    # drop leading 14-digit timestamp
    return re.sub(r"^\d{14}", "", name)


def group_by_suffix(path):
    out = {}
    for f in os.listdir(path):
        if f.endswith(".pickle"):
            out[f] = f
    return out


def compare_pickle_collections(legacy_dir, new_dir, label):
    leg = sorted(os.listdir(legacy_dir))
    new = sorted(os.listdir(new_dir))
    problems = []
    if leg != new:
        problems.append(f"{label}: file set differs only-legacy={sorted(set(leg)-set(new))} only-new={sorted(set(new)-set(leg))}")
    for f in sorted(set(leg) & set(new)):
        import pickle
        try:
            lp = pickle.load(open(legacy_dir / f, "rb"))
            npd = pickle.load(open(new_dir / f, "rb"))
            assert_pickle_equivalent(npd, lp, f"pickle.{label}.{f}")
        except (Mismatch, Exception) as e:
            problems.append(f"{label}/{f}: {e}")
    return problems


def compare_global_dir(legacy_dir, new_dir, label):
    # group by stripped suffix, compare the most recent in each dir
    def latest(path):
        out = {}
        for f in os.listdir(path):
            out.setdefault(strip_ts(f), []).append(f)
        return {k: max(v) for k, v in out.items()}

    leg = latest(legacy_dir)
    new = latest(new_dir)
    problems = []
    if set(leg) != set(new):
        problems.append(f"{label}: suffix set differs only-legacy={sorted(set(leg)-set(new))} only-new={sorted(set(new)-set(leg))}")
    for suf in sorted(set(leg) & set(new)):
        lf = leg[suf]
        nf = new[suf]
        if suf.endswith(".json"):
            import json
            a = json.load(open(new_dir / nf))
            b = json.load(open(legacy_dir / lf))
            if a != b:
                problems.append(f"{label}/{suf}: JSON differs")
        elif suf.endswith(".npy"):
            a = np.load(new_dir / nf)
            b = np.load(legacy_dir / lf)
            if not np.array_equal(a, b):
                problems.append(f"{label}/{suf}: NPY differs")
    return problems


def fig_sets(legacy_dir, new_dir, label):
    def figs(path):
        return {strip_ts(f): f for f in os.listdir(path) if f.endswith(".png")}
    leg = figs(legacy_dir)
    new = figs(new_dir)
    problems = []
    if set(leg) != set(new):
        problems.append(f"{label}: fig set differs only-legacy={sorted(set(leg)-set(new))} only-new={sorted(set(new)-set(leg))}")
    if not (set(leg) == set(new)):
        return problems
    # compare dimensions for matching figure kinds
    for key in sorted(set(leg) & set(new)):
        try:
            import matplotlib.image as mpimg
            a = mpimg.imread(legacy_dir / leg[key]).shape
            b = mpimg.imread(new_dir / new[key]).shape
            if a != b:
                problems.append(f"{label}/{key}: dimensions {a} != {b}")
        except Exception as e:
            problems.append(f"{label}/{key}: {e}")
    return problems


def main():
    problems = []
    # Level 1: single-fit pickles
    problems += compare_pickle_collections(ORACLE / "TB_fit_result", NEW / "single_process" / "TB_fit_result", "TB single")
    problems += compare_pickle_collections(ORACLE / "EC_fit_result", NEW / "single_process" / "EC_fit_result", "EC single")
    # Level 2/3: global outputs
    problems += compare_global_dir(ORACLE / "tb_logs", NEW / "tb_logs", "tb_logs")
    problems += compare_global_dir(ORACLE / "ec_logs", NEW / "ec_logs", "ec_logs")
    # figures
    problems += fig_sets(ORACLE / "single_fit_fig", NEW / "single_process" / "single_fit_fig", "single fig")
    problems += fig_sets(ORACLE / "tb_logs", NEW / "tb_logs", "tb fig")
    problems += fig_sets(ORACLE / "ec_logs", NEW / "ec_logs", "ec fig")
    if problems:
        print(f"FAIL ({len(problems)}):")
        for p in problems:
            print("  -", p)
        return 1
    print("PASS: all 09 legacy vs new outputs equivalent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
