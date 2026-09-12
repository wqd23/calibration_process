# -*- coding:utf-8 -*-
"""Unit coverage for the GRID-N1 workflow module (selection + TB enumeration)."""

from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest
import yaml

from calibration_process import pipeline
from calibration_process.workflows.registry import get_workflow
from calibration_process.workflows.versions import vN1


def test_n1_registered_for_all_subversions():
    for ver in vN1.VERSIONS:
        assert get_workflow(ver) is vN1


def test_selection_only_for_neutron(monkeypatch):
    assert vN1.selection("N1-Gamma-Am241", "tb") is None
    assert vN1.selection("N1-Neutron", "tb") is None

    monkeypatch.setattr(vN1, "NEUTRON_CCM_MIN", 0.5)
    selkey, fn = vN1.selection("N1-Neutron", "neutron")
    assert selkey == vN1.NEUTRON_SELKEY
    frames = {"data_max": np.array([100, 200, 300]),
              "data_base": np.array([0, 0, 0]),
              "data_ccm": np.array([65535, 0, 32768])}
    mask, extra = fn(frames)
    assert list(mask) == [True, False, True]
    assert extra["ccm"][0] == 1.0
    assert extra["amp"][1] == 200.0


def _rt():
    return SimpleNamespace(payload=SimpleNamespace(tb=SimpleNamespace(science_dir="raw_data")))


def _make_raw(tmp_path):
    raw = tmp_path / "raw_data"
    raw.mkdir()
    for name in ("0C-265-290-139.event.dat", "m20C-265-125.event.dat",
                 "0C-285-CI-137.event.dat", "0To10C-265-165.event.dat",
                 "not_event.dat"):
        (raw / name).write_bytes(b"")


def test_enumerate_tb(tmp_path):
    _make_raw(tmp_path)
    gagg = vN1.enumerate_measurements("N1-Gamma-Am241", "tb", _rt(), tmp_path)
    ids = [r["id"] for r in gagg]
    # 0C scan -> 8 standard biases, m20C single -> 1; CI/To excluded
    assert len(ids) == 9
    assert "GAGG_0C_265" in ids and "GAGG_0C_290" in ids and "GAGG_m20C_265" in ids
    seg = {r["id"]: r["metadata"]["seg_bias"] for r in gagg}
    assert seg["GAGG_0C_265"] == 265
    assert seg["GAGG_m20C_265"] is None

    clyc = vN1.enumerate_measurements("N1-Gamma-Na22", "tb", _rt(), tmp_path)
    clyc_ids = [r["id"] for r in clyc]
    # 0C scan gains the extra 281 segment; m20C single stays
    assert "CLYC_0C_281" in clyc_ids and "CLYC_m20C_265" in clyc_ids


def test_enumerate_matches_committed_fit_range():
    for ver in ("N1-Gamma-Am241", "N1-Gamma-Na22"):
        if not (Path("data") / ver / "raw_data").exists():
            pytest.skip("N1 raw_data not linked")
        rt = pipeline.load_rt(ver)
        ids = {r["id"] for r in vN1.enumerate_measurements(ver, "tb", rt, Path("data") / ver)}
        fr = yaml.safe_load(
            (Path("src/calibration_process/configs") / ver / "fit_range_tb.yaml").read_text()
        )["measurements"]
        assert ids == set(fr), f"{ver}: {sorted(ids ^ set(fr))[:6]}"
