# -*- coding:utf-8 -*-
"""End-to-end run of the explicit new workflow (09, 12B).

Runs the full formal workflow (single fit + global fit for TB and EC) into a
temporary output root, then compares every output to the frozen legacy oracle
(``.oracle/<ver>``).  This exercises the pipeline, the shared stages and the
version workflow, and provides the Level 1/2/3 regression evidence.

Skipped when the raw data or the frozen oracle are absent.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from compare import (  # noqa: E402
    assert_json_equivalent,
    assert_npy_equivalent,
    assert_pickle_equivalent,
)

SRC = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC / "calibration_process"))
from calibration_process.cli import main as cli_main  # noqa: E402

VERSIONS = ["09", "12B"]
# expected single-fit pickle counts per version
COUNTS = {"09": (48, 16), "12B": (54, 18)}


def test_all_versions_match_oracle(tmp_path):
    for VER in VERSIONS:
        oracle = Path(".oracle") / VER
        if not (oracle / "TB_fit_result").exists():
            pytest.skip(f".oracle/{VER} not present")
        out = tmp_path / VER
        assert cli_main(["all", VER, "-o", str(out)]) == 0

        tb_n, ec_n = COUNTS[VER]
        assert len(list((out / "single_process/TB_fit_result").glob("*.pickle"))) == tb_n
        assert len(list((out / "single_process/EC_fit_result").glob("*.pickle"))) == ec_n

        import pickle
        for sub in ("TB_fit_result", "EC_fit_result"):
            leg = oracle / sub
            new = out / "single_process" / sub
            for f in sorted(os_listdir(leg)):
                assert_pickle_equivalent(
                    pickle.load(open(new / f, "rb")), pickle.load(open(leg / f, "rb")),
                    f"pickle.{sub}.{f}",
                )
        _compare_global(out / "tb_logs", oracle / "tb_logs")
        _compare_global(out / "ec_logs", oracle / "ec_logs")
        _compare_figs(out / "single_process/single_fit_fig", oracle / "single_fit_fig")
        _compare_figs(out / "tb_logs", oracle / "tb_logs")
        _compare_figs(out / "ec_logs", oracle / "ec_logs")


def os_listdir(p):
    import os
    return sorted(os.listdir(p))


def _compare_global(new_dir, leg_dir):
    import os
    import re

    def latest(path):
        out = {}
        for f in os.listdir(path):
            out[re.sub(r"^\d{14}", "", f)] = f
        return out

    leg, new = latest(leg_dir), latest(new_dir)
    assert set(leg) == set(new), f"{leg_dir}: {set(leg) ^ set(new)}"
    import json
    for suf in sorted(set(leg) & set(new)):
        lf, nf = leg[suf], new[suf]
        if suf.endswith(".json"):
            assert_json_equivalent(json.load(open(new_dir / nf)),
                                   json.load(open(leg_dir / lf)), "json." + suf)
        elif suf.endswith(".npy"):
            assert_npy_equivalent(np.load(new_dir / nf), np.load(leg_dir / lf), suf)


def _compare_figs(new_dir, leg_dir):
    import os
    import re
    import matplotlib.image as mpimg

    def figs(path):
        return {re.sub(r"^\d{14}", "", f): f for f in os.listdir(path) if f.endswith(".png")}

    leg, new = figs(leg_dir), figs(new_dir)
    assert set(leg) == set(new), f"{leg_dir}: {set(leg) ^ set(new)}"
    for key in sorted(set(leg) & set(new)):
        assert mpimg.imread(leg_dir / leg[key]).shape == mpimg.imread(new_dir / new[key]).shape
