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


def _stub_rt(**ec_kw):
    """A minimal runtime stub with the fields single_run_spec needs."""
    from types import SimpleNamespace
    ec = SimpleNamespace(
        reader=ec_kw.get("reader", "xray"),
        xray_reader=ec_kw.get("xray_reader"),
        xray_single_file=ec_kw.get("xray_single_file", False),
        xray_config_file=ec_kw.get("xray_config_file"),
        time_cut=ec_kw.get("time_cut"),
        bin_width=ec_kw.get("bin_width", 10),
        adc_max=ec_kw.get("adc_max", 16384.0),
        channel_count=ec_kw.get("channel_count", 4),
    )
    return SimpleNamespace(
        data_dir=Path("data/05B"),
        payload=SimpleNamespace(ec=ec),
        corr=[lambda t, b: 1.0] * 4,
        fit_ranges={"ec_xray": {"15keV": [[0, 10], [0, 10], [0, 10], [0, 10]]}},
        bkg_forms={"ec_xray": {"15keV": "lin"}},
        fit_range=lambda branch, key: {"15keV": [[0, 10], [0, 10], [0, 10], [0, 10]]}[key],
        bkg_form=lambda branch, key: "lin",
    )


def test_single_run_spec_xray_single_file_05b():
    from calibration_process.config_schema import ManifestEntry
    rt = _stub_rt(xray_single_file=True, xray_reader="xray",
                  xray_config_file="data/05B/raw_data/x/config.json",
                  time_cut={"15keV.dat": [[1, 2], [3, 4], [5, 6], [7, 8]]})
    m = ManifestEntry(id="15keV", branch="ec_xray",
                      science_files=["raw_data/x/15keV.dat"])
    spec = stages.single_run_spec(rt, "ec_xray", m)
    assert spec.read_config.ending == "xray"
    assert spec.read_config.config_file.endswith("config.json")
    # bkg time cut is the cyclically-rotated per-channel range
    assert spec.bkg_read_config.time_cut == [[3, 4], [5, 6], [7, 8], [1, 2]]
    assert spec.read_config.time_cut == [[1, 2], [3, 4], [5, 6], [7, 8]]
    assert spec.fit_config.bkg_form == "lin"


def test_time_cut_helpers():
    from types import SimpleNamespace
    pb = SimpleNamespace(time_cut={"a.dat": [[1, 2], [3, 4], [5, 6], [7, 8]]})
    assert stages._time_cut(pb, "a.dat") == [[1, 2], [3, 4], [5, 6], [7, 8]]
    assert stages._time_cut(pb, "b.dat") is None
    assert stages._bkg_time_cut(pb, "a.dat") == [[3, 4], [5, 6], [7, 8], [1, 2]]
    assert stages._bkg_time_cut(pb, "b.dat") is None
    assert stages._time_cut(SimpleNamespace(time_cut=None), "a.dat") is None


def test_pad_to4_and_pad_group_4ch():
    assert stages._pad_to4([1, 2, 3]) == [1, 2, 3, 1]
    assert stages._pad_to4([1, 2, 3, 4]) == [1, 2, 3, 4]
    groups = [[{"b": 1}], [{"b": 1}, {"b": 2}, {"b": 3}], [{"b": 1}, {"b": 2}, {"b": 3}, {"b": 4}]]
    padded = stages._pad_group_4ch(groups)
    assert len(padded[0]) == 4
    assert len(padded[1]) == 4
    assert len(padded[2]) == 4


def test_resolution_fit_methods():
    import numpy as np
    e = np.array([10.0, 20.0, 30.0, 40.0])
    r = np.array([0.20, 0.10, 0.07, 0.05])
    re_ = np.array([0.02, 0.01, 0.007, 0.005])
    popt, perr = stages._resolution_fit("exprfit", e, r, re_)
    assert len(popt) == 3
    popt, perr = stages._resolution_fit("lmfit", e, r, re_)
    assert len(popt) == 3
    with pytest.raises(ValueError):
        stages._resolution_fit("bogus", e, r, re_)


def test_build_ec_points_respects_channel_count():
    from calibration_process.config_schema import ManifestEntry
    from types import SimpleNamespace
    rt = SimpleNamespace(
        payload=SimpleNamespace(ec=SimpleNamespace(channel_count=3)),
        energies={"20": 20.0},
    )
    m = ManifestEntry(id="20", branch="ec_xray")
    fp = _fake_fp([
        {"b": 1, "b_err": 0.1, "a": 1, "a_err": 0.1, "c": 1, "c_err": 0.1,
         "resolution": 0.1, "resolution_err": 0.01, "rate": 1, "rate_err": 0.1,
         "redchi": 1, "ndf": 5, "success": True, "qa_flag": "ok"},
        {"b": 2, "b_err": 0.1, "a": 1, "a_err": 0.1, "c": 1, "c_err": 0.1,
         "resolution": 0.1, "resolution_err": 0.01, "rate": 1, "rate_err": 0.1,
         "redchi": 1, "ndf": 5, "success": True, "qa_flag": "ok"},
        {"b": 3, "b_err": 0.1, "a": 1, "a_err": 0.1, "c": 1, "c_err": 0.1,
         "resolution": 0.1, "resolution_err": 0.01, "rate": 1, "rate_err": 0.1,
         "redchi": 1, "ndf": 5, "success": True, "qa_flag": "ok"},
        {"b": 4, "b_err": 0.1, "a": 1, "a_err": 0.1, "c": 1, "c_err": 0.1,
         "resolution": 0.1, "resolution_err": 0.01, "rate": 1, "rate_err": 0.1,
         "redchi": 1, "ndf": 5, "success": True, "qa_flag": "ok"},
    ])
    per_ch = stages.build_ec_points(rt, [(m, fp)], "xray")
    assert len(per_ch) == 3                 # channel_count caps to 3
    for ch_list in per_ch:
        assert len(ch_list) == 1            # only channel ch contributed
