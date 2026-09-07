# -*- coding:utf-8 -*-
"""End-to-end run of the explicit new workflow for version 09.

Runs the full formal workflow (single fit + global fit for TB and EC) into a
temporary output root, then compares every output to the frozen legacy oracle
(``.oracle/09``).  This exercises the pipeline, the shared stages and the
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

VER = "09"
ORACLE = Path(".oracle") / VER
pytestmark = pytest.mark.skipif(
    not (ORACLE / "TB_fit_result").exists(),
    reason=".oracle/09 not present (run the legacy pipeline first)",
)


def test_full_pipeline_matches_oracle(tmp_path):
    out = tmp_path / "out"
    assert cli_main(["all", VER, "-o", str(out)]) == 0

    # file set / counts
    tb_pickles = out / "single_process/TB_fit_result"
    ec_pickles = out / "single_process/EC_fit_result"
    assert len(list(tb_pickles.glob("*.pickle"))) == 48
    assert len(list(ec_pickles.glob("*.pickle"))) == 16

    # Level 1: pickles
    import pickle
    for sub in ("TB_fit_result", "EC_fit_result"):
        leg = ORACLE / sub
        new = out / "single_process" / sub
        for f in sorted(os_listdir(leg)):
            assert_pickle_equivalent(
                pickle.load(open(new / f, "rb")), pickle.load(open(leg / f, "rb")),
                f"pickle.{sub}.{f}",
            )

    # Level 3: global outputs
    _compare_global(out / "tb_logs", ORACLE / "tb_logs")
    _compare_global(out / "ec_logs", ORACLE / "ec_logs")

    # figure sets + dimensions
    _compare_figs(out / "single_process/single_fit_fig", ORACLE / "single_fit_fig")
    _compare_figs(out / "tb_logs", ORACLE / "tb_logs")
    _compare_figs(out / "ec_logs", ORACLE / "ec_logs")


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
