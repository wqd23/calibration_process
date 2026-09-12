# -*- coding:utf-8 -*-
"""Unit coverage for the GRID-N1 workflow module (selection + enumeration)."""

import json

import numpy as np

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


def test_enumerate_reads_file_map(tmp_path):
    class RT:
        config_root = tmp_path

    assert vN1.enumerate_measurements("N1-Neutron", "neutron", RT(), tmp_path) == []
    (tmp_path / "neutron_file_map.json").write_text(json.dumps([
        {"id": "m1", "science_files": ["raw_data/a.event.dat"],
         "hk_files": ["raw_data/a.hk"], "metadata": {"mode": "wf"}}]))
    out = vN1.enumerate_measurements("N1-Neutron", "neutron", RT(), tmp_path)
    assert len(out) == 1
    assert out[0]["id"] == "m1"
    assert out[0]["branch"] == "neutron"
    assert out[0]["science_files"] == ["raw_data/a.event.dat"]
    assert out[0]["metadata"] == {"mode": "wf"}
