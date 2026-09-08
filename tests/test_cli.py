# -*- coding:utf-8 -*-
"""CLI subcommand coverage (light, read-mostly)."""

import sys
from pathlib import Path

SRC = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC / "calibration_process"))
from calibration_process.cli import main, build_parser  # noqa: E402

VER = "09"
ID = "0825_0C_265_4m_0x00C5.txt"


def test_list(capsys):
    assert main(["list", VER, "tb"]) == 0
    out = capsys.readouterr().out
    assert ID in out


def test_config_show(capsys):
    assert main(["config", VER, "tb", ID]) == 0
    out = capsys.readouterr().out
    assert "fit" in out.lower()
    assert ID in out


def test_fit_one_to_tmp(tmp_path):
    out = tmp_path / "out"
    assert main(["fit-one", VER, "tb", ID, "-o", str(out)]) == 0
    assert (out / "single_process/TB_fit_result/0825_0C_265_4m_0x00C5.pickle").exists()
    assert (out / "single_process/single_fit_fig/0825_0C_265_4m_0x00C5.png").exists()


def test_parser_has_new_commands():
    p = build_parser()
    choices = p._subparsers._group_actions[0].choices
    for c in ("discover", "list", "fit-one", "fit", "global", "all", "config"):
        assert c in choices


def test_all_version_nocache_output_root(tmp_path):
    # light: only verify output_root envelope wiring (no full run)
    from calibration_process import pipeline
    o = pipeline.output_root(VER, str(tmp_path / "x"))
    assert o == tmp_path / "x"


def test_check_cli(capsys):
    assert main(["check", "09"]) == 0
    assert "READY" in capsys.readouterr().out


def test_config_not_found_raises():
    import pytest as pt
    with pt.raises(SystemExit):
        main(["config", "09", "tb", "zzz_not_here"])


def test_global_unknown_branch_raises():
    import pytest as pt
    with pt.raises(SystemExit):
        main(["global", "09", "nope"])
