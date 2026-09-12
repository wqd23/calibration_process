# -*- coding:utf-8 -*-
"""Unit coverage for scripts/merge_tb.py."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import merge_tb as mt  # noqa: E402


def _coeff(tag):
    return {"G0": tag, "k": 1.0, "V0": 24.0, "b": -35.0, "c": -1000.0}


def _write(tmp_path, name, slots):
    p = tmp_path / name
    p.write_text(json.dumps(slots))
    return str(p)


def test_merge_by_channel(tmp_path):
    a = _write(tmp_path, "a.json", [None, _coeff("a1"), _coeff("a2"), None])
    b = _write(tmp_path, "b.json", [_coeff("b0"), None, None, _coeff("b3")])
    merged = mt.merge_tb([(a, [1, 2]), (b, [0, 3])])
    assert merged[0]["G0"] == "b0"
    assert merged[1]["G0"] == "a1"
    assert merged[2]["G0"] == "a2"
    assert merged[3]["G0"] == "b3"


def test_merge_partial_leaves_null(tmp_path):
    a = _write(tmp_path, "a.json", [None, _coeff("a1"), None, None])
    merged = mt.merge_tb([(a, [1])])
    assert merged[1]["G0"] == "a1"
    assert merged[0] is None and merged[2] is None and merged[3] is None


def test_merge_duplicate_channel_raises(tmp_path):
    a = _write(tmp_path, "a.json", [_coeff("a0")] + [None] * 3)
    b = _write(tmp_path, "b.json", [_coeff("b0")] + [None] * 3)
    with pytest.raises(ValueError):
        mt.merge_tb([(a, [0]), (b, [0])])


def test_parse_source_and_cli(tmp_path):
    assert mt._parse_source("x.json:1,2") == ("x.json", [1, 2])
    with pytest.raises(ValueError):
        mt._parse_source("x.json")

    a = _write(tmp_path, "a.json", [None, _coeff("a1"), _coeff("a2"), None])
    b = _write(tmp_path, "b.json", [_coeff("b0"), None, None, _coeff("b3")])
    out = tmp_path / "merged.json"
    rc = mt.main(["-o", str(out), "--from", f"{a}:1,2", "--from", f"{b}:0,3"])
    assert rc == 0
    merged = json.loads(out.read_text())
    assert [v["G0"] for v in merged] == ["b0", "a1", "a2", "b3"]
