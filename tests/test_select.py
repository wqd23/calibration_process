# -*- coding:utf-8 -*-
"""Unit coverage for the L1->L2 event-selection hook and its sidecar cache."""

import numpy as np
import pytest

from lib_reader import EventTable
from lib_reader import l1_cache


def test_event_table_row():
    frames = {"a": np.array([1, 2, 3]), "wf": np.arange(6).reshape(3, 2)}
    t = EventTable(frames)
    assert len(t) == 3
    row = t.row(1)
    assert row["a"] == 2
    assert list(row["wf"]) == [2, 3]


def _masked_amps(result):
    return list(np.asarray(result["amp"]))


def test_apply_selection_mask_and_sidecar(tmp_path, monkeypatch):
    monkeypatch.setattr(l1_cache, "_cache_root", lambda ver: tmp_path)
    frames = {"amp": np.array([1, 2, 3, 4]), "ch": np.array([0, 1, 0, 1])}
    calls = {"n": 0}

    def fn(f):
        calls["n"] += 1
        return f["amp"] > 2, {"score": f["amp"] * 10.0}

    masked = l1_cache.apply_selection("T", "r", "/x.dat", {}, "sel@1", fn, frames)
    assert calls["n"] == 1
    assert _masked_amps(masked) == [3, 4]
    assert list(masked["ch"]) == [0, 1]
    assert list(tmp_path.rglob("select__*.parquet"))

    # a repeat with the same selkey reuses the sidecar (fn is not called again)
    masked2 = l1_cache.apply_selection("T", "r", "/x.dat", {}, "sel@1", fn, frames)
    assert calls["n"] == 1
    assert _masked_amps(masked2) == [3, 4]

    # a different selkey is a different selection and recomputes
    l1_cache.apply_selection("T", "r", "/x.dat", {}, "sel@2", fn, frames)
    assert calls["n"] == 2


def test_apply_selection_identity_and_bad_mask(tmp_path, monkeypatch):
    monkeypatch.setattr(l1_cache, "_cache_root", lambda ver: tmp_path)
    frames = {"a": np.array([1, 2, 3])}
    masked = l1_cache.apply_selection(
        "T", "r", "/y.dat", {}, "all", lambda f: np.ones(len(f["a"]), bool), frames
    )
    assert list(masked["a"]) == [1, 2, 3]
    with pytest.raises(ValueError):
        l1_cache.apply_selection(
            "T", "r", "/z.dat", {}, "bad", lambda f: np.array([True]), frames
        )
