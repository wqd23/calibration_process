# -*- coding:utf-8 -*-
"""Differential regression: manifest set + config resolution vs legacy (09, 12B).

These tests are fast (no data reading) and prove the strongest equivalence
argument: the explicit workflow resolves the *same* read/bkg/spectrum/fit
configs as the historical ``file_config`` for every selected measurement.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
SRC = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC / "calibration_process"))
from calibration_process import manifest as man  # noqa: E402
from calibration_process.pipeline import load_rt, _manifest_path  # noqa: E402
from calibration_process.workflows import common as stages  # noqa: E402
from legacy_ops import legacy_tb, legacy_ec, legacy_file_config  # noqa: E402

VERSIONS = ["03B", "04", "05B", "07", "10B", "11B", "09", "12B"]
BRANCHES = ["tb", "ec_source", "ec_xray"]


@pytest.mark.parametrize("ver", VERSIONS)
@pytest.mark.parametrize("branch", BRANCHES)
def test_manifest_set_matches_legacy(ver, branch):
    manifest = man.load_manifest(_manifest_path(ver, branch))
    ids = [m.id for m in manifest.measurements]
    if branch == "tb":
        legacy_ids = legacy_tb(ver).files
    elif branch == "ec_source":
        legacy_ids = legacy_ec(ver).src_list
    else:
        legacy_ids = list(legacy_ec(ver).x_list)
    assert set(ids) == set(legacy_ids), f"{ver}/{branch}: {sorted(set(ids) ^ set(legacy_ids))}"
    assert len(ids) == len(set(ids)), "duplicated measurement id"


@pytest.mark.parametrize("ver", VERSIONS)
@pytest.mark.parametrize("branch", BRANCHES)
def test_config_resolution_matches_legacy(ver, branch):
    rt = load_rt(ver)
    manifest = man.load_manifest(_manifest_path(ver, branch))
    for m in manifest.measurements:
        legacy_cfg = legacy_file_config(ver, branch, m.id)
        new_spec = stages.single_run_spec(rt, branch, m)
        _assert_spec_equivalent(legacy_cfg, new_spec, f"{ver}/{branch}/{m.id}")


def _assert_spec_equivalent(legacy_cfg, new_spec, label):
    from dataclasses import fields as dc_fields

    lr, lb, ls, lf = legacy_cfg
    nr, nb, ns, nf = (new_spec.read_config, new_spec.bkg_read_config,
                      new_spec.spectrum_config, new_spec.fit_config)

    def d(c):
        return {f.name: getattr(c, f.name) for f in dc_fields(c)}

    for name, a, b in (("read", lr, nr), ("bkg", lb, nb)):
        if isinstance(a, list):
            assert len(a) == len(b), label
            for i, (x, y) in enumerate(zip(a, b)):
                assert d(x) == d(y), f"{label} {name}[{i}]"
        else:
            assert d(a) == d(b), f"{label} {name}"
    for k in ("adc_max", "bin_width", "rate_style", "corr_style"):
        assert d(ls)[k] == d(ns)[k], f"{label} spectrum.{k}"
    assert len(d(ls)["corr"]) == len(d(ns)["corr"]), f"{label} corr len"
    assert d(lf) == d(nf), f"{label} fit"
