# -*- coding:utf-8 -*-
"""Self-contained golden regression for the unified B/C readers.

Runs the *new* unified readers on small committed real raw samples
(``tests/golden/reader/<sample>/raw/``) and asserts the ``(sci, tel)`` output
equals the frozen *legacy* output (``expected.npz`` + ``structure.json``).  No
``raw_data`` / ``.oracle`` dependency.  Covers the distinct reader paths: C hex
(07/04), B waveform noUdp (05B normal), B waveform UDP with little-endian HK
(05B xray), B waveform UDP (03B src) and B feature UDP with the per-file time cut
(03B xray).  Regenerate with ``scripts/gen_reader_golden.py``.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from lib_reader.reader07.frame_adapter import single_read04, single_read07
from lib_reader.reader05.frame_adapter import (
    single_read05b_normal,
    single_read05b_xray,
    single_read03b,
    src_read03b,
)
from lib_reader.readerN1.read import single_readN1

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "golden" / "reader"

# sample -> (reader, rundata filename, reader args after path)
CASES = [
    ("07_hex", single_read07, "bnu_Am241_10cm_15min_220423180307_COM3-Data.txt", ()),
    ("04_hex", single_read04, "210501124621_COM7_tb_-20C_27p0V_4m_5cm-Data.txt", ()),
    ("05B_normal", single_read05b_normal, "TB_0C_275_rundata2022-04-11-02-29-33.dat", ()),
    ("05B_xray", single_read05b_xray, "XM_100_rundata2022-04-25-13-46-36.dat", ("",)),
    ("03B_src", src_read03b, "src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat", ("",)),
    ("03B_xray", single_read03b, "jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat", ("",)),
]


def _reconstruct(flat, structure):
    out = {}
    for section in ("sci", "tel"):
        d = {}
        for k, kind in structure[section].items():
            if kind == "list4":
                d[k] = [flat[f"{section}.{k}.{i}"] for i in range(4)]
            elif kind == "empty":
                d[k] = []
            else:
                d[k] = flat[f"{section}.{k}"]
        out[section] = d
    return out["sci"], out["tel"]


def _assert_pair_equal(a, b, tag):
    assert type(a) is type(b), (tag, type(a), type(b))
    for section, idx in (("sci", 0), ("tel", 1)):
        assert set(a[idx].keys()) == set(b[idx].keys()), (tag, section)
        for k in a[idx]:
            x, y = a[idx][k], b[idx][k]
            if isinstance(x, list) or isinstance(y, list):
                assert isinstance(x, list) and isinstance(y, list), (tag, section, k)
                assert len(x) == len(y), (tag, section, k)
                for i, (u, v) in enumerate(zip(x, y)):
                    u, v = np.asarray(u), np.asarray(v)
                    assert u.dtype == v.dtype and np.array_equal(u, v), (tag, section, k, i)
            else:
                x, y = np.asarray(x), np.asarray(y)
                assert x.dtype == y.dtype and np.array_equal(x, y), (tag, section, k)


@pytest.mark.parametrize("sample,reader,rundata,args", CASES)
def test_reader_golden(sample, reader, rundata, args):
    raw_dir = OUT / sample / "raw"
    if not (raw_dir / rundata).exists():
        pytest.skip(f"golden raw sample missing: {raw_dir / rundata}")

    structure = json.loads((OUT / sample / "structure.json").read_text())
    with np.load(OUT / sample / "expected.npz") as npz:
        expected = _reconstruct({k: npz[k] for k in npz.files}, structure)

    produced = reader(str(raw_dir / rundata), *args)
    _assert_pair_equal(produced, expected, f"reader.{sample}")


N1_EVENT = "0degC-28.5-191.event.dat"
N1_HK = "0degC-28.5-ecu_113.hk"


def test_reader_n1_golden():
    # GRID-N1 has no legacy implementation; this golden is self-consistent
    raw_dir = OUT / "n1_wf" / "raw"
    if not (raw_dir / N1_EVENT).exists():
        pytest.skip(f"golden raw sample missing: {raw_dir / N1_EVENT}")
    structure = json.loads((OUT / "n1_wf" / "structure.json").read_text())
    with np.load(OUT / "n1_wf" / "expected.npz") as npz:
        expected = _reconstruct({k: npz[k] for k in npz.files}, structure)
    produced = single_readN1(str(raw_dir / N1_EVENT), hk_path=str(raw_dir / N1_HK), mode="wf")
    _assert_pair_equal(produced, expected, "reader.n1_wf")
