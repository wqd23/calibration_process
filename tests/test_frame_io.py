# -*- coding:utf-8 -*-
"""Tests for the byte-source layer (binary / hex-text)."""
from pathlib import Path

import pytest

from lib_reader.frame_io import load_binary, load_hex_text

ROOT = Path(__file__).resolve().parents[1]

D07 = ROOT / "data/07/raw_data/北师大上胶后补标定/bnu_Am241_10cm_15min_220423180307_COM3-Data.txt"
D05 = ROOT / "data/05B/raw_data/test/TB_0C_275_rundata2022-04-11-02-29-33.dat"


def test_load_binary(tmp_path):
    p = tmp_path / "x.dat"
    p.write_bytes(bytes([0, 1, 127, 255]))
    assert load_binary(p).tolist() == [0, 1, 127, 255]


def test_load_hex_text(tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("AA bb 02\n0F 00")
    assert load_hex_text(p).tolist() == [0xAA, 0xBB, 0x02, 0x0F, 0x00]


def test_hex_text_real_head():
    if not D07.exists():
        pytest.skip(f"raw data missing: {D07}")
    assert load_hex_text(D07)[:4].tolist() == [0xAA, 0xBB, 0xCC, 0x02]


def test_load_binary_real():
    if not D05.exists():
        pytest.skip(f"raw data missing: {D05}")
    b = load_binary(D05)
    assert b.dtype.name == "uint8" and b.shape[0] > 0
