# -*- coding:utf-8 -*-
"""Unit coverage for the shared stage helpers (no data reading)."""

import sys
from pathlib import Path
import numpy as np
import pytest

SRC = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC / "calibration_process"))
from calibration_process.workflows import common as stages  # noqa: E402
from calibration_process.config_schema import ManifestEntry  # noqa: E402


def _fake_fp(fit_results):
    class FP:
        pass
    fp = FP()
    fp.fit_result = fit_results
    fp.tel = {"tempSipm": [np.array([20.0, 22.0])] * 4,
              "bias": [np.array([28.0, 28.2])] * 4}
    return fp


def test_to_single_fit_result_with_none_channel():
    m = ManifestEntry(id="x", branch="tb")
    fr = [{
        "a": 1.0, "a_err": 0.1, "b": 100.0, "b_err": 0.5, "c": 2.0, "c_err": 0.2,
        "resolution": 0.05, "resolution_err": 0.005, "rate": 10.0, "rate_err": 1.0,
        "redchi": 1.2, "ndf": 50, "success": True, "qa_flag": "ok",
        "boundary_hit": [],
    }, None]
    fp = _fake_fp(fr)
    out = stages.to_single_fit_result(m, fp)
    assert out[0].peak_center == 100.0
    assert out[0].channel == 0
    assert out[1] is None


def test_rotate_bkg_fixed_and_unknown():
    reads = list(range(4))
    assert stages._rotate_bkg(reads, "circle", 4) == [1, 2, 3, 0]
    assert stages._rotate_bkg(reads, "fixed", 4) == [1, 2, 0, 0]
    with pytest.raises(ValueError):
        stages._rotate_bkg(reads, "bogus", 4)


def test_single_run_spec_unknown_branch():
    from calibration_process.config_schema import ManifestEntry
    from types import SimpleNamespace
    m = ManifestEntry(id="x", branch="tb")
    rt = SimpleNamespace(data_dir=Path("data/09"))
    with pytest.raises(ValueError):
        stages.single_run_spec(rt, "nope", m)


def test_apply_bias_filter():
    data = np.array([[1, 1, 1, 1, 27.0, 1], [1, 1, 1, 1, 28.0, 1]])
    assert stages._apply_bias_filter(data, None).shape[0] == 2
    assert stages._apply_bias_filter(data, 27.5).shape[0] == 1


def test_build_tb_points_skips_none_fit():
    rt = None
    m = ManifestEntry(id="m", branch="tb")
    fp = _fake_fp([None, {
        "b": 100.0, "b_err": 0.5, "a": 1, "a_err": 0.1, "c": 2, "c_err": 0.2,
        "resolution": 0.05, "resolution_err": 0.005, "rate": 1, "rate_err": 0.1,
        "redchi": 1, "ndf": 5, "success": True, "qa_flag": "ok"}])
    pts = stages.build_tb_points(rt, [(m, fp)])
    assert len(pts[0]) == 0
    assert len(pts[1]) == 1
    assert pts[1][0].measurement_id == "m"


def test_build_ec_points_channel_cap():
    from calibration_process.config_schema import ManifestEntry
    # channel_count ignored via runtime here; just verify fields
    rt = None
    m = ManifestEntry(id="15keV", branch="ec_xray")
    fp = _fake_fp([{"b": 10.0, "b_err": 0.1, "a": 1, "a_err": 0.1, "c": 2, "c_err": 0.2,
                    "resolution": 0.1, "resolution_err": 0.01, "rate": 1, "rate_err": 0.1,
                    "redchi": 1, "ndf": 5, "success": True, "qa_flag": "ok"}])
    # call with a stub runtime holding energies + channel_count
    class Rt:
        payload = type("P", (), {"ec": type("E", (), {"channel_count": 4})()})()
        energies = {"15keV": 15.0}
    rt = Rt()
    pts = stages.build_ec_points(rt, [(m, fp)], "xray")
    assert len(pts[0]) == 1
    assert pts[0][0].source_kind == "xray"
