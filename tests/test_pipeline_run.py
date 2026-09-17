# -*- coding:utf-8 -*-
"""End-to-end run of the explicit new workflow (09, 12B).

Runs the full formal workflow (single fit + global fit for TB and EC) into a
temporary output root, then compares every output to the frozen legacy oracle
(``<oracle>/<ver>``).  This exercises the pipeline, the shared stages and the
version workflow, and provides the Level 1/2/3 regression evidence.

The oracle root defaults to ``.oracle`` in the repository (usually a symlink
to the shared fixture directory) and can be overridden with the
``CALIB_ORACLE_DIR`` environment variable.  Skipped when the raw data or the
frozen oracle are absent.
"""

import os
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

# Versions whose reader now reports the SiPM-side bias (no 499-ohm
# series-resistor drop): this also shifts the HK stable-epoch selection, so the
# whole telemetry differs from the frozen legacy oracle.  These versions still
# run the pipeline and the output counts are checked, but the byte-level oracle
# diff is skipped until the hardware convention is confirmed (see
# docs/intermediate_data.md section 6.1).  09 is untouched and still strict.
PENDING_BIAS = {"12B"}

ORACLE_ROOT = Path(os.environ.get("CALIB_ORACLE_DIR", ".oracle"))


def test_all_versions_match_oracle(tmp_path):
    for VER in VERSIONS:
        oracle = ORACLE_ROOT / VER
        if not (oracle / "TB_fit_result").exists():
            pytest.skip(f"{oracle} not present")
        out = tmp_path / VER
        assert cli_main(["all", VER, "-o", str(out)]) == 0

        tb_n, ec_n = COUNTS[VER]
        assert len(list((out / "single_process/TB_fit_result").glob("*.pickle"))) == tb_n
        assert len(list((out / "single_process/EC_fit_result").glob("*.pickle"))) == ec_n

        if VER in PENDING_BIAS:
            print(f"NOTE: skip {VER} legacy oracle diff: SiPM-side bias "
                  f"convention pending (docs/intermediate_data.md 6.1)")
            continue

        import pickle
        for sub in ("TB_fit_result", "EC_fit_result"):
            leg = oracle / sub
            new = out / "single_process" / sub
            for f in sorted(os_listdir(leg)):
                if not f.endswith(".pickle"):
                    # the oracle may also carry L3 artefacts (fit.json /
                    # spectrum.parquet); only the pickles are compared here
                    continue
                try:
                    legacy_obj = pickle.load(open(leg / f, "rb"))
                except (ModuleNotFoundError, ImportError, AttributeError) as e:
                    # Frozen legacy EC pickles are dill-bound to the removed
                    # ``operation`` module and cannot be unpickled in the new
                    # package.  EC scientific identity is instead validated by
                    # the global ec_logs json/npy comparison below, so the
                    # per-measurement pickle is skipped with a recorded reason.
                    print(f"NOTE: skip legacy {sub}/{f} (pickle coupling: {e})")
                    continue
                assert_pickle_equivalent(
                    pickle.load(open(new / f, "rb")), legacy_obj,
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
    # ``new`` figures must all have a legacy counterpart (and match its shape).
    # The archived oracle may contain stale figures from older dev runs (e.g.
    # a dropped 12B "15" point), so legacy-only leftovers are tolerated.
    assert set(new) <= set(leg), f"{leg_dir}: new-only {sorted(set(new) - set(leg))}"
    for key in sorted(set(leg) & set(new)):
        assert mpimg.imread(leg_dir / leg[key]).shape == mpimg.imread(new_dir / new[key]).shape
