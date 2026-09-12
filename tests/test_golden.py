# -*- coding:utf-8 -*-
"""Small, self-contained golden regression tests for the explicit workflow.

Two families, both free of any ``raw_data`` / ``.oracle`` dependency at test
time (they run on a fresh clone with only the committed configs and the small
``data/*/single_process/*.json`` references):

A. Orchestration behaviours of ``single_run_spec`` for the distinct
   version branches (tb / ec_source / ec_xray single-file / multi-file with
   circle vs fixed background rotation, hk passthrough).

B. Global-fit stages (``global_tb`` / ``global_ec``) reconstructed from
   committed, reduced real points (``tests/golden/<ver>/points_*.json``) and
   matched back to the committed golden coefficients (which are the frozen
   legacy oracle).  Each representative version exercises a distinct behaviour:

      09  tb curvefit + ec polyfit(4ch)
      11B tb lmfit + ec channel_count 3
      12B tb bias_min_filter
      03B ec lmfit resolution

   ``04`` (exprfit resolution) and ``10B`` (channel_count 3) are NOT
   golden-tested here because the legacy run produced no EC global output for
   them (no authoritative reference); they are unit-covered in test_stages.py.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from compare import assert_json_equivalent, assert_npy_equivalent  # noqa: E402

SRC = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC / "calibration_process"))

from calibration_process import manifest as man  # noqa: E402
from calibration_process.products import ECPoint, TBPoint  # noqa: E402
from calibration_process.workflows import common as stages  # noqa: E402

DATA = Path("data")
CONFIG = Path("src/calibration_process/configs")
GOLDEN = Path(__file__).parent / "golden"

TB_GOLDEN = ["09", "11B", "12B", "GRIDN1/GAGG", "GRIDN1/CLYC"]
EC_GOLDEN = ["09", "03B", "11B", "GRIDN1/EC"]


def _rt(ver):
    from calibration_process.pipeline import load_rt
    return load_rt(ver)


def _man(ver, branch):
    return man.load_manifest(CONFIG / ver / f"{branch}_manifest.yaml")


def _first(ver, branch):
    return man.filtered_measurements(_man(ver, branch))[0]


# --------------------------------------------------------------------------- #
# A. single_run_spec orchestration behaviours
# --------------------------------------------------------------------------- #
def _spec(ver, branch, measurement=None):
    from calibration_process.workflows import common as stages
    rt = _rt(ver)
    m = measurement or _first(ver, branch)
    return rt, m, stages.single_run_spec(rt, branch, m)


def _abspath(ver, rel):
    return str(DATA / ver / rel)


def test_spec_tb_resolution():
    rt, m, spec = _spec("09", "tb")
    assert spec.read_config.path == _abspath("09", m.science_files[-1])
    assert spec.read_config.ending == rt.payload.tb.reader
    assert spec.spectrum_config.bin_width == rt.payload.tb.bin_width
    assert spec.spectrum_config.adc_max == rt.payload.tb.adc_max
    assert spec.fit_config.fit_range == rt.fit_range("tb", m.id)
    assert spec.fit_config.bkg_form == rt.bkg_form("tb", m.id)


def test_spec_tb_hk_passthrough():
    # 12B TB measurement that carries an hk file + sci_half/hk_bias metadata
    rt, m, spec = _spec("12B", "tb", measurement=_man("12B", "tb").measurements[1])
    assert spec.read_config.kwarg["hk_path"] == _abspath("12B", m.hk_files[-1])
    assert spec.read_config.kwarg["sci_half"] == "first"
    assert spec.read_config.kwarg["hk_bias"] == 27.0


def test_spec_ec_source_bkg_and_corr():
    rt, m, spec = _spec("09", "ec_source")
    assert spec.read_config.path == _abspath("09", m.science_files[-1])
    assert spec.bkg_read_config.path == _abspath("09", m.aux_files[-1])
    assert len(spec.spectrum_config.corr) == rt.payload.ec.channel_count


def test_spec_ec_xray_single_file():
    # 05B: one multi-channel file, science and background share the path
    rt, m, spec = _spec("05B", "ec_xray")
    assert spec.read_config.path == _abspath("05B", m.science_files[-1])
    assert spec.bkg_read_config.path == _abspath("05B", m.science_files[-1])
    assert spec.read_config.ending == (rt.payload.ec.xray_reader or rt.payload.ec.reader)


def test_spec_ec_xray_multi_circle_rotation():
    # 09: ch i background comes from channel (i+1)%4
    rt, m, spec = _spec("09", "ec_xray")
    reads, bkg = spec.read_config, spec.bkg_read_config
    assert len(reads) == 4 and len(bkg) == 4
    for i in range(4):
        assert bkg[i].path == reads[(i + 1) % 4].path
    for i, f in enumerate(m.science_files):
        assert reads[i].path == _abspath("09", f)


def test_spec_ec_xray_multi_fixed_rotation():
    # 12B: legacy fixed rotation [ch1, ch2, ch0, ch0]
    rt, m, spec = _spec("12B", "ec_xray")
    reads, bkg = spec.read_config, spec.bkg_read_config
    assert len(reads) == 4 and len(bkg) == 4
    for i, src in enumerate([1, 2, 0, 0]):
        assert bkg[i].path == reads[src].path


def test_spec_all_versions_configured():
    # every version resolvable for all its branches (smoke over the matrix)
    for ver in ["03B", "04", "05B", "07", "09", "10B", "11B", "12B"]:
        for branch, fname in [("tb", "tb_manifest"), ("ec_source", "ec_source_manifest"),
                              ("ec_xray", "ec_xray_manifest")]:
            if (CONFIG / ver / f"{fname}.yaml").exists():
                _spec(ver, branch)


# --------------------------------------------------------------------------- #
# B. global-fit golden (reconstructed points -> coefficients)
# --------------------------------------------------------------------------- #
def _load_points_tb(ver):
    data = json.load(open(GOLDEN / ver / "points_tb.json"))
    per_channel = [[] for _ in range(4)]
    for p in data:
        per_channel[p["channel"]].append(TBPoint(**p))
    return per_channel


def _load_points_ec(ver):
    data = json.load(open(GOLDEN / ver / "points_ec.json"))
    src, xray = [[] for _ in range(4)], [[] for _ in range(4)]
    for p in data:
        kind = p.pop("source_kind")
        (src if kind == "src" else xray)[p["channel"]].append(ECPoint(source_kind=kind, **p))
    return src, xray


def _latest(path, suffix):
    cands = [f for f in path.iterdir() if f.name.endswith(suffix)]
    assert cands, f"no {suffix} in {path}"
    return max(cands)


@pytest.mark.parametrize("ver", TB_GOLDEN)
def test_global_tb_golden(tmp_path, ver):
    rt = _rt(ver)
    per_channel = _load_points_tb(ver)
    stages.global_tb(rt, per_channel, tmp_path / "tb_logs")
    produced = json.load(open(_latest(tmp_path / "tb_logs", "temp_bias_fit.json")))
    golden = json.load(open(GOLDEN / ver / "tb_coeff.json"))
    assert_json_equivalent(produced, golden, f"tb.{ver}")


def test_global_tb_neutron_selfconsistent(tmp_path):
    """N1 neutron runs: fixed bias, so the 5-parameter TB fit is degenerate.

    This is a *self-consistent* anchor (not a legacy/physics baseline): it locks
    the current fit output so later refactors are compared bit-for-bit.  The
    frozen points come from the 260322 neutron temperature scans at 28.5 V.
    """
    rt = _rt("GRIDN1/Neutron")
    per_channel = _load_points_tb("GRIDN1/Neutron")
    stages.global_tb(rt, per_channel, tmp_path / "tb_logs")
    produced = json.load(open(_latest(tmp_path / "tb_logs", "temp_bias_fit.json")))
    golden = json.load(open(GOLDEN / "GRIDN1/Neutron" / "tb_coeff.json"))
    assert_json_equivalent(produced, golden, "tb.GRIDN1/Neutron")


@pytest.mark.parametrize("ver", EC_GOLDEN)
def test_global_ec_golden(tmp_path, ver):
    rt = _rt(ver)
    src, xray = _load_points_ec(ver)
    stages.global_ec(rt, src, xray, tmp_path / "ec_logs")
    nch = rt.payload.ec.channel_count
    for ch in range(nch):
        produced = json.load(open(_latest(tmp_path / "ec_logs", f"ec_coef_sci_ch{ch}.json")))
        golden = json.load(open(GOLDEN / ver / f"ec_coeff_ch{ch}.json"))
        assert_json_equivalent(produced, golden, f"ec.{ver}.ch{ch}")
        prod_npy = np.load(_latest(tmp_path / "ec_logs", f"ec_data_ch{ch}.npy"))
        gold_npy = np.load(GOLDEN / ver / f"ec_data_ch{ch}.npy")
        assert_npy_equivalent(prod_npy, gold_npy, f"ec.{ver}.ch{ch}")


def test_golden_points_file_is_small():
    # sanity: committed golden must stay tiny (guard against bloat).  The B/C
    # reader golden adds small truncated raw samples + frozen npz, so the budget
    # is larger than the points-only golden (~240 KB).
    total = sum(f.stat().st_size for f in GOLDEN.rglob("*") if f.is_file())
    assert total < 3 * 1024 * 1024, f"golden grew to {total} bytes"
