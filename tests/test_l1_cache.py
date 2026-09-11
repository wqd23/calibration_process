# -*- coding:utf-8 -*-
"""Tests for the L1 parquet cache (reader11 / reader12 raw-parse layer)."""
import json
import os
import importlib
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from addict import Dict

from lib_reader import l1_cache
from lib_reader.l1_cache import (
    deserialize,
    read_frame,
    serialize,
    write_frame,
)

ROOT = Path(__file__).resolve().parents[1]

MODS = {
    "reader12": "lib_reader.reader12.read",
    "reader11": "lib_reader.reader11.read",
    "reader10": "lib_reader.reader10.read",
}

D12 = ROOT / "data/12B/raw_data/温度偏压/074/TB_m0C_270.dat"
H12 = D12.with_suffix(".hk")
D11 = (
    ROOT
    / "data/11B/raw_data/tempbias/温度偏压0～10摄氏度/251_0_Cs137_28.5_observe.dat"
)
H11 = D11.with_name("251_0_Cs137_28.5.hk")
D10 = ROOT / "data/10B/raw_data/src_data/073_observe_Cs137.dat"
H10 = D10.with_name("073_hk_Cs137.dat")

# (ver, reader, dat, hk, mode, rmod); mode=None => readSci takes no mode arg
CASES = [
    pytest.param("12B", "12b", D12, H12, "ft", "reader12", id="reader12"),
    pytest.param("11B", "11b", D11, H11, "wf", "reader11", id="reader11"),
    pytest.param("10B", "10b", D10, H10, None, "reader10", id="reader10"),
]


def _require(path: Path):
    if not path.exists():
        pytest.skip(f"raw data missing: {path}")
    return path


def _sci_post(raw: Dict) -> dict:
    """A reduced slice of single_readXX's sci post-processing: amp + channel split."""
    amp = np.asarray(raw["data_max"]) - np.asarray(raw["data_base"]) / 4.0
    raw["amp"] = amp
    raw["timestampEvt"] = raw["timestamp"]
    n = raw["data_max"].shape[0]
    out = {}
    for k, v in raw.items():
        arr = np.asarray(v)
        if arr.shape[0] == n and k != "channel_n":
            out[k] = [arr[raw["channel_n"] == i] for i in range(4)]
        else:
            out[k] = arr
    return out


def _assert_dicts_equal(a: Dict, b: Dict):
    assert type(a) is type(b), (type(a), type(b))
    assert set(a.keys()) == set(b.keys())
    for k in a:
        av, bv = np.asarray(a[k]), np.asarray(b[k])
        assert np.asarray(av).dtype == np.asarray(bv).dtype, (k, av.dtype, bv.dtype)
        assert np.array_equal(av, bv), k


def _module(rmod):
    return importlib.import_module(MODS[rmod])


def _sci_read(m, dat, mode, **kw):
    if mode is None:
        return m.readSci(dat, **kw)
    return m.readSci(dat, mode=mode, **kw)


def _sci_impl(m, dat, mode):
    f = m._readSci_impl
    if mode is None:
        return f(dat)
    return f(dat, mode=mode)


def _sci_record(ver, raw_ptr):
    index = json.loads((l1_cache._cache_root(ver) / "index.json").read_text())
    return next(
        r for r in index if r["kind"] == "sci" and r["raw_ptr"] == os.fspath(raw_ptr)
    )


@pytest.mark.parametrize("ver,reader,dat,hk,mode,rmod", CASES)
def test_serializer_roundtrip(ver, reader, dat, hk, mode, rmod, tmp_path):
    _require(dat)
    impl = getattr(_module(rmod), "_readSci_impl")
    raw = impl(dat) if mode is None else impl(dat, mode=mode)
    cols, dtype, rows, shapes = serialize(raw)
    p = tmp_path / "rt.parquet"
    write_frame(cols, p)
    back = read_frame(p, {"dtypes": dtype, "shapes": shapes})
    _assert_dicts_equal(Dict(back), raw)
    assert rows == len(np.asarray(raw["data_max"]))


@pytest.mark.parametrize("ver,reader,dat,hk,mode,rmod", CASES)
def test_l1_reconstruct_equals_direct(ver, reader, dat, hk, mode, rmod):
    _require(dat)
    m = _module(rmod)
    sci_raw = _sci_read(m, dat, mode, overwrite_cache=True)
    direct = _sci_post(Dict({k: np.asarray(v) for k, v in sci_raw.items()}))

    rec = _sci_record(ver, dat)
    cols = read_frame(ROOT / rec["events"])
    recovered = _sci_post(deserialize(cols))

    for k in direct:
        d, r = direct[k], recovered[k]
        if isinstance(d, list):
            assert len(d) == len(r)
            for da, ra in zip(d, r):
                assert np.asarray(da).dtype == np.asarray(ra).dtype
                assert np.array_equal(np.asarray(da), np.asarray(ra))
        else:
            assert np.asarray(d).dtype == np.asarray(r).dtype
            assert np.array_equal(np.asarray(d), np.asarray(r))


def test_l1_miss_then_hit(monkeypatch, tmp_path):
    m = _module("reader12")
    calls = {"n": 0}

    def fake_impl(path, mode="ft"):
        calls["n"] += 1
        return Dict(
            {
                "data_max": np.array([10, 20, 30], dtype=np.uint32),
                "data_base": np.array([4, 4, 4], dtype=np.uint32),
                "channel_n": np.array([0, 1, 2], dtype=np.uint32),
                "crc_check": np.array([True, False, True]),
            }
        )

    monkeypatch.setattr(m, "_readSci_impl", fake_impl)
    fake_path = str(tmp_path / "fake_sci.dat")

    m.readSci(fake_path, mode="ft")
    assert calls["n"] == 1
    m.readSci(fake_path, mode="ft")
    assert calls["n"] == 1, "second call must hit the L1 cache"

    index = json.loads((l1_cache._cache_root("12B") / "index.json").read_text())
    rec = next(r for r in index if r["raw_ptr"] == fake_path and r["kind"] == "sci")
    assert (l1_cache._cache_root("12B") / "cache" / rec["key"] / "meta.json").exists()
    assert (l1_cache._cache_root("12B") / "cache" / rec["key"] / "events.parquet").exists()


def test_overwrite_cache_forces_reparse(monkeypatch, tmp_path):
    m = _module("reader12")
    calls = {"n": 0}

    def fake_impl(path, mode="ft"):
        calls["n"] += 1
        return Dict(
            {
                "data_max": np.array([1], dtype=np.uint32),
                "data_base": np.array([0], dtype=np.uint32),
                "channel_n": np.array([0], dtype=np.uint32),
            }
        )

    monkeypatch.setattr(m, "_readSci_impl", fake_impl)
    fake_path = str(tmp_path / "fake_overwrite.dat")
    m.readSci(fake_path, mode="ft", overwrite_cache=True)
    m.readSci(fake_path, mode="ft", overwrite_cache=True)
    assert calls["n"] == 2


def test_corrupt_falls_back():
    m = _module("reader12")
    real_path = _require(D12)

    out = m.readSci(real_path, mode="ft")
    assert isinstance(out, Dict)

    rec = _sci_record("12B", real_path)
    events = ROOT / rec["events"]
    events.write_bytes(b"\x00garbage-not-a-parquet")
    out2 = m.readSci(real_path, mode="ft")
    assert isinstance(out2, Dict)
    assert np.array_equal(out2["data_max"], out["data_max"])


def test_index_external_read():
    """A pure polars+json consumer (no lib_reader import) reads the sci parquet."""
    m = _module("reader12")
    m.readSci(D12, mode="ft", overwrite_cache=True)
    rec = _sci_record("12B", D12)
    df = pl.read_parquet(ROOT / rec["events"])
    amp = df["data_max"].to_numpy() - df["data_base"].to_numpy() / 4.0
    assert amp.shape[0] > 0
