# -*- coding:utf-8 -*-
"""Tests for the processed-output L1 parquet cache (four-channel stacking).

Covers the monolithic readers (04/07/09/03B/05B) whose final ``(sci, tel)``
output has jagged list-of-4 channel quantities.  Missing raw data -> skip.
"""
import json
from pathlib import Path

import numpy as np
import pytest

import lib_reader
from lib_reader import l1_cache
from lib_reader.l1_cache import (
    deserialize_processed,
    serialize_processed,
    with_l1_cache_processed,
)

ROOT = Path(__file__).resolve().parents[1]

D04 = ROOT / "data/04/raw_data/20210504_source_GRID04/210504162829_COM6_src_Co60_10m_10cm-Data.txt"
D05N = ROOT / "data/05B/raw_data/test/TB_0C_275_rundata2022-04-11-02-29-33.dat"
D05X = ROOT / "data/05B/raw_data/X光机实验-天格/XM_100_rundata2022-04-25-13-46-36.dat"
D03S = ROOT / "data/03B/raw_data/20210504_source_03B/src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat"
D03X = ROOT / "data/03B/raw_data/20210429_Xray_03B/jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat"
D07 = ROOT / "data/07/raw_data/北师大上胶后补标定/bnu_Am241_10cm_15min_220423180307_COM3-Data.txt"
D09 = ROOT / "data/09/raw_data/Xray/15keV_ch0_0x00EF_2m.txt"

# (id, l1 func, args)
CASES = [
    pytest.param("04", "single_read04", (D04,), id="04"),
    pytest.param("05B_normal", "single_read05b_normal", (D05N,), id="05B_normal"),
    pytest.param("05B_xray", "single_read05b_xray", (D05X, ""), id="05B_xray"),
    pytest.param("03B_src", "src_read03b", (D03S, ""), id="03B_src"),
    pytest.param("03B_xray", "single_read03b", (D03X, ""), id="03B_xray"),
    pytest.param("07", "single_read07", (D07,), id="07"),
    pytest.param("09", "single_read09", (D09,), id="09"),
]


def _require(case):
    for p in case:
        if isinstance(p, Path) and not p.exists():
            pytest.skip(f"raw data missing: {p}")


def _call(fn, args):
    return fn(*[str(a) if isinstance(a, Path) else a for a in args])


def _assert_equal(a, b, path=""):
    if isinstance(a, list) or isinstance(b, list):
        assert isinstance(a, list) and isinstance(b, list), path
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _assert_equal(x, y, f"{path}[{i}]")
    elif isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        assert isinstance(a, np.ndarray) and isinstance(b, np.ndarray), path
        assert a.dtype == b.dtype, (path, a.dtype, b.dtype)
        assert np.array_equal(a, b), path
    else:
        assert a == b, path


def _assert_pair_equal(a, b, tag):
    assert type(a) is type(b), (tag, type(a), type(b))
    for section, idx in (("sci", 0), ("tel", 1)):
        assert set(a[idx].keys()) == set(b[idx].keys()), (tag, section)
        for k in a[idx]:
            _assert_equal(a[idx][k], b[idx][k], f"{tag}.{section}.{k}")


@pytest.mark.parametrize("ver,l1fn,args", CASES)
def test_roundtrip_serializer(ver, l1fn, args):
    _require(args)
    result = _call(getattr(lib_reader, l1fn), args)
    section_d, section_meta = serialize_processed(result)
    sci, tel = deserialize_processed(section_d, section_meta)
    _assert_pair_equal((sci, tel), result, f"roundtrip.{ver}")


def _fake_result():
    sci = {
        "amp": [
            np.arange(3, dtype=np.uint32),
            np.arange(5, dtype=np.uint32),
            np.array([], dtype=np.uint32),
            np.arange(2, dtype=np.uint32),
        ],
        "timestampEvt": [
            np.arange(3, dtype=np.float64),
            np.arange(5, dtype=np.float64),
            np.array([], dtype=np.float64),
            np.arange(2, dtype=np.float64),
        ],
        "sciNum": [np.array([], dtype=np.int64) for _ in range(4)],
        "effectiveCount": np.arange(4, dtype=np.int64),
        "timeCorrect": [],
    }
    tel = {
        "tempSipm": [np.arange(2, dtype=np.float64) for _ in range(4)],
        "timestamp": np.arange(2, dtype=np.float64),
        "telNum": [],
    }
    return sci, tel


def test_miss_then_hit(tmp_path):
    calls = {"n": 0}

    @with_l1_cache_processed(ver="L1TEST", reader="fake", kind="single")
    def read(path):
        calls["n"] += 1
        return _fake_result()

    p = str(tmp_path / "fake.dat")
    a = read(p)
    assert calls["n"] == 1
    b = read(p)
    assert calls["n"] == 1, "second call must hit the cache"
    _assert_pair_equal(a, b, "miss_hit")

    index = json.loads((l1_cache._cache_root("L1TEST") / "index.json").read_text())
    rec = next(r for r in index if r["raw_ptr"] == p and r["kind"] == "single")
    assert (l1_cache._cache_root("L1TEST") / "cache" / rec["key"] / "meta.json").exists()


def test_overwrite_forces_recompute(tmp_path):
    calls = {"n": 0}

    @with_l1_cache_processed(ver="L1TEST", reader="fake", kind="single")
    def read(path):
        calls["n"] += 1
        return _fake_result()

    p = str(tmp_path / "over.dat")
    read(p, overwrite_cache=True)
    read(p, overwrite_cache=True)
    assert calls["n"] == 2


def test_corrupt_falls_back(tmp_path):
    calls = {"n": 0}

    @with_l1_cache_processed(ver="L1TEST", reader="fake", kind="single")
    def read(path):
        calls["n"] += 1
        return _fake_result()

    p = str(tmp_path / "corrupt.dat")
    a = read(p)
    index = json.loads((l1_cache._cache_root("L1TEST") / "index.json").read_text())
    rec = next(r for r in index if r["raw_ptr"] == p and r["kind"] == "single")
    files = list((l1_cache._cache_root("L1TEST") / "cache" / rec["key"]).glob("*.parquet"))
    assert files
    files[0].write_bytes(b"\x00garbage")
    b = read(p)
    assert calls["n"] == 2
    _assert_pair_equal(a, b, "corrupt")
