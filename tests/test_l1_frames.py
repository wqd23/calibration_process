# -*- coding:utf-8 -*-
"""L1 faithful-frame smoke test on the committed reader golden samples.

This is a structural lock (not a legacy baseline): it checks that the unified
``lib_reader.read_frames`` returns the expected L1 columns — notably a boolean
``crc_check`` — with consistent row counts, without adding new fixtures.  The
scientific (L2) equivalence to legacy stays covered by ``test_reader_golden``.
"""
from pathlib import Path

import numpy as np
import pytest

import lib_reader

OUT = Path(__file__).resolve().parent / "golden" / "reader"

N1_EVENT = "0degC-28.5-191.event.dat"
N1_HK = "0degC-28.5-ecu_113.hk"

# sample -> (ver, kind, raw file name, read_frames kwargs)
CASES = [
    ("04_hex", "04", "sci", "210501124621_COM7_tb_-20C_27p0V_4m_5cm-Data.txt", {}),
    ("07_hex", "07", "sci", "bnu_Am241_10cm_15min_220423180307_COM3-Data.txt", {}),
    ("05B_normal", "05B", "sci", "TB_0C_275_rundata2022-04-11-02-29-33.dat", {}),
    ("05B_normal", "05B", "hk", "TB_0C_275_HK2022-04-11-02-29-33.dat", {"ending": "normal"}),
    ("05B_xray", "05B", "sci", "XM_100_rundata2022-04-25-13-46-36.dat", {}),
    ("05B_xray", "05B", "hk", "XM_100_HK2022-04-25-13-46-36.dat", {"ending": "x_ray"}),
    ("03B_src", "03B", "sci", "src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat",
     {"feature_mode": False}),
    ("03B_src", "03B", "hk", "src_Cs137_12m_10cm_HK2021-05-05-12-12-27.dat", {}),
    ("03B_src", "03B", "tl", "src_Cs137_12m_10cm_TimeLine2021-05-05-12-12-27.dat", {}),
    ("03B_xray", "03B", "sci", "jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat",
     {"feature_mode": True}),
    ("n1_wf", "N1-Neutron", "sci", N1_EVENT, {"mode": "wf"}),
    ("n1_wf", "N1-Neutron", "hk", N1_HK, {}),
]


@pytest.mark.parametrize("sample,ver,kind,fname,kwargs", CASES)
def test_read_frames_structure(sample, ver, kind, fname, kwargs):
    path = OUT / sample / "raw" / fname
    if not path.exists():
        pytest.skip(f"golden raw sample missing: {path}")
    frames = lib_reader.read_frames(str(path), ver, kind, **kwargs)
    assert isinstance(frames, dict) and frames, (sample, kind)
    rows = None
    if kind == "sci":
        # every science frame carries the parse-time CRC verdict
        assert "crc_check" in frames, (sample, kind)
        crc = np.asarray(frames["crc_check"])
        assert crc.dtype == np.bool_, (sample, kind, crc.dtype)
        rows = len(crc)
        assert rows > 0
    for key, value in frames.items():
        arr = np.asarray(value)
        assert arr.shape[0] > 0, (sample, kind, key)
        if rows is not None:
            assert arr.shape[0] == rows, (sample, kind, key)
