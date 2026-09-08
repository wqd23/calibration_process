# -*- coding:utf-8 -*-
"""Unit coverage for the new config / manifest / orchestration layer."""

import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))
SRC = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC / "calibration_process"))
from calibration_process import pipeline  # noqa: E402
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


def test_registry_unknown_version_raises():
    with pytest.raises(KeyError):
        get_workflow("99X")


@pytest.mark.parametrize("ver", ["03B", "04", "05B", "07", "10B", "11B", "09", "12B"])
@pytest.mark.parametrize("branch", ["tb", "ec_source", "ec_xray"])
def test_versions_enumerate_matches_manifest(ver, branch):
    # the committed manifest is the frozen, legacy-validated measurement set;
    # discover->enumerate must reproduce it exactly
    rt = pipeline.load_rt(ver)
    data_dir = Path("data") / ver
    wf = get_workflow(ver)
    got = {r["id"] for r in wf.enumerate_measurements(ver, branch, rt, data_dir)}
    expected = {m.id for m in load_manifest(
        Path("src/calibration_process/configs") / ver / f"{branch}_manifest.yaml").measurements}
    assert got == expected, f"{ver}/{branch}: {sorted(got ^ expected)}"


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
