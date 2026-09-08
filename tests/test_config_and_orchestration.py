# -*- coding:utf-8 -*-
"""Unit coverage for the new config / migration / orchestration layer."""

import sys
from pathlib import Path

import pytest
import yaml

SRC = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC / "calibration_process"))
from calibration_process import migration, pipeline  # noqa: E402
from calibration_process import manifest as man  # noqa: E402
from calibration_process.manifest import discover, load_manifest, filtered_measurements  # noqa: E402
from calibration_process.workflows.registry import get_workflow  # noqa: E402
from calibration_process.runtime import load_runtime  # noqa: E402
from calibration_process.config_schema import (  # noqa: E402
    PayloadSchema,
    AnalysisSchema,
    FitRangeSet,
    ManifestSchema,
)

VER = "09"
CONFIG_ROOT = Path("src/calibration_process/configs") / VER


def test_payload_analysis_schemas_load():
    p = PayloadSchema.model_validate(yaml.safe_load(open(CONFIG_ROOT / "payload.yaml")))
    assert p.version == "09"
    assert p.ec.channel_count == 4
    a = AnalysisSchema.model_validate(yaml.safe_load(open(CONFIG_ROOT / "analysis.yaml")))
    assert a.ec_xray.default_bkg is None
    assert a.ec_source.background_overrides["0827_10C_285_Co60_20m_0x010F.txt"] == "gaus"


def test_fit_range_schema_strict():
    FitRangeSet.model_validate({"measurements": {"x": [[0, 1]] * 4}})
    # a per-channel None range is allowed (12B uses un-fitted channels)
    FitRangeSet.model_validate({"measurements": {"x": [[0, 1], None, [0, 1], None]}})
    with pytest.raises(Exception):
        FitRangeSet.model_validate({"measurements": {"x": [[0, 1], [0, 1]]}})
    with pytest.raises(Exception):
        FitRangeSet.model_validate({"measurements": {"x": [[1, 0]] * 4}})
    with pytest.raises(Exception):
        FitRangeSet.model_validate({"measurements": {"x": [[0, 1]] * 4}, "extra": 1})


def test_manifest_schema_duplicate_id_rejected():
    with pytest.raises(Exception):
        ManifestSchema.model_validate({
            "version": "09", "branch": "tb",
            "measurements": [
                {"id": "a", "branch": "tb"},
                {"id": "a", "branch": "tb"},
            ],
        })


def test_fit_range_schema_accepts_real():
    d = yaml.safe_load(open(CONFIG_ROOT / "fit_range_tb.yaml"))
    fr = FitRangeSet.model_validate(d)
    assert len(fr.measurements) == 48


def test_manifest_discover_and_load():
    for branch in ("tb", "ec_source", "ec_xray"):
        rt = pipeline.load_rt(VER)
        m = discover(VER, branch, rt, Path("data") / VER)
        assert len(m.measurements) > 0
        assert m.version == VER
        # round trip through dump/load
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "m.yaml"
            man.dump_manifest(m, p)
            m2 = load_manifest(p)
            assert len(m2.measurements) == len(m.measurements)
        assert len(filtered_measurements(m)) == len(m.measurements)


def test_versions_v09_enumerate_matches_legacy():
    from calibration_process import operation as op
    rt = pipeline.load_rt(VER)
    data_dir = Path("data") / VER
    wf = get_workflow(VER)
    tb_ids = [r["id"] for r in wf.enumerate_measurements(VER, "tb", rt, data_dir)]
    leg_tb = op.TB_operation_09(
        path="data/09/raw_data/temp_bias", fit_range="data/09/single_process/fit_range.json",
        save_path="x", save_fig_path="x", result_path="x")
    assert set(tb_ids) == set(leg_tb.files)
    ec_op = op.EC_operation_09(
        tb_result_path="data/09/single_process/20260301164328_temp_bias_fit.json",
        fit_range="data/09/single_process/fit_range.json",
        energy="data/09/single_process/ec_energy.json",
        bkg_form="data/09/single_process/bkg_form.json",
        x_path="data/09/raw_data/Xray", src_path="data/09/raw_data/src",
        save_path="x", save_fig_path="x", result_path="x")
    src_ids = [r["id"] for r in wf.enumerate_measurements(VER, "ec_source", rt, data_dir)]
    assert set(src_ids) == set(ec_op.src_list)
    x_ids = [r["id"] for r in wf.enumerate_measurements(VER, "ec_xray", rt, data_dir)]
    assert set(x_ids) == set(ec_op.x_list)


def test_registry_unknown_version_raises():
    with pytest.raises(KeyError):
        get_workflow("99X")


def test_versions_v12b_enumerate_matches_legacy():
    from calibration_process import operation as op
    rt = pipeline.load_rt("12B")
    data_dir = Path("data") / "12B"
    wf = get_workflow("12B")
    tb_ids = [r["id"] for r in wf.enumerate_measurements("12B", "tb", rt, data_dir)]
    leg_tb = op.TB_operation_12B(
        path="data/12B/raw_data", fit_range="data/12B/single_process/fit_range.json",
        save_path="x", save_fig_path="x", result_path="x",
        file_map="data/12B/single_process/tb_file_map.json")
    assert set(tb_ids) == set(leg_tb.files)
    ec_op = op.EC_operation_12B(
        tb_result_path="data/12B/single_process/20260824134441_temp_bias_fit.json",
        fit_range="data/12B/single_process/fit_range.json",
        energy="data/12B/single_process/ec_energy.json",
        bkg_form="data/12B/single_process/bkg_form.json",
        x_path="data/12B/raw_data/X光机/072", src_path="data/12B/raw_data/放射源/072",
        save_path="x", save_fig_path="x", result_path="x")
    src_ids = [r["id"] for r in wf.enumerate_measurements("12B", "ec_source", rt, data_dir)]
    assert set(src_ids) == set(ec_op.src_list)
    x_ids = [r["id"] for r in wf.enumerate_measurements("12B", "ec_xray", rt, data_dir)]
    assert set(x_ids) == set(ec_op.x_list)


def test_migration_writes_valid_yaml_12b(tmp_path):
    root = tmp_path / "configs" / "12B"
    migration.migrate_version("12B", config_root=root)
    p = PayloadSchema.model_validate(yaml.safe_load(open(root / "payload.yaml")))
    assert p.tb.bias_min_filter == 27.25
    assert p.tb.tb_file_map == "single_process/tb_file_map.json"
    assert p.ec.xray_bkg_rotation == "fixed"
    assert p.ec.energy_map["0611_Co60_25min_240f0032.dat"] == 1332.492
    tb_fr = FitRangeSet.model_validate(yaml.safe_load(open(root / "fit_range_tb.yaml")))
    assert len(tb_fr.measurements) == 54
    x_fr = FitRangeSet.model_validate(yaml.safe_load(open(root / "fit_range_ec_xray.yaml")))
    assert len(x_fr.measurements) == 14


def test_migration_writes_valid_yaml(tmp_path):
    root = tmp_path / "configs" / VER
    migration.migrate_version(VER, config_root=root)
    p = PayloadSchema.model_validate(yaml.safe_load(open(root / "payload.yaml")))
    assert p.ec.energy_map["0827_10C_285_Co60_20m_0x010F.txt"] == 1332.0
    assert len(FitRangeSet.model_validate(
        yaml.safe_load(open(root / "fit_range_tb.yaml"))).measurements) == 48
    assert len(FitRangeSet.model_validate(
        yaml.safe_load(open(root / "fit_range_ec_xray.yaml"))).measurements) == 13


def test_runtime_corr_from_tb_ref():
    rt = load_runtime(VER, CONFIG_ROOT, Path("data") / VER)
    assert len(rt.corr) == 4
    for i in range(4):
        f = rt.corr[i]
        assert abs(f(25.0, 28.5) - 1.0) < 1e-9  # corr is 1 at the reference point


def test_channel_use_semantics():
    from calibration_process.workflows.common import channel_use
    m = load_manifest(CONFIG_ROOT / "tb_manifest.yaml").measurements[0]
    assert channel_use(m, 0) is True
    # measurement-level disable
    from calibration_process.config_schema import ManifestEntry
    disabled = ManifestEntry(id="x", branch="tb", use=False, channels=None)
    assert channel_use(disabled, 0) is False
    ch_disabled = ManifestEntry(id="y", branch="tb", channels={"2": {"use": False}})
    assert channel_use(ch_disabled, 2) is False
    assert channel_use(ch_disabled, 0) is True
