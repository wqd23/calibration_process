# -*- coding:utf-8 -*-
"""Factories that build the legacy operation objects for a version.

Used by the differential regression tests to compare the new explicit config
resolution against the historical ``file_config``-produced configs.  These must
mirror the constructor arguments from ``operation.py`` for each version.
"""

from pathlib import Path

from calibration_process import operation as op

DATA = Path("data")


def _sp(ver):
    return DATA / ver / "single_process"


def legacy_tb(ver):
    if ver == "09":
        return op.TB_operation_09(
            path=str(DATA / ver / "raw_data/temp_bias"),
            fit_range=str(_sp(ver) / "fit_range.json"),
            save_path=str(_sp(ver) / "TB_fit_result"),
            save_fig_path=str(_sp(ver) / "single_fit_fig"),
            result_path=str(DATA / ver / "tb_logs"),
        )
    if ver == "12B":
        return op.TB_operation_12B(
            path=str(DATA / ver / "raw_data"),
            fit_range=str(_sp(ver) / "fit_range.json"),
            save_path=str(_sp(ver) / "TB_fit_result"),
            save_fig_path=str(_sp(ver) / "single_fit_fig"),
            result_path=str(DATA / ver / "tb_logs"),
            file_map=str(_sp(ver) / "tb_file_map.json"),
        )
    raise ValueError(f"no legacy TB op factory for {ver}")


def legacy_ec(ver):
    sp = _sp(ver)
    if ver == "09":
        tb_ref = "20260301164328_temp_bias_fit.json"
        x_path = "raw_data/Xray"
        src_path = "raw_data/src"
    elif ver == "12B":
        tb_ref = "20260824134441_temp_bias_fit.json"
        x_path = "raw_data/X光机/072"
        src_path = "raw_data/放射源/072"
    else:
        raise ValueError(f"no legacy EC op factory for {ver}")
    return op.EC_operation_05B if False else _make_ec(ver, sp, tb_ref, x_path, src_path)


def _make_ec(ver, sp, tb_ref, x_path, src_path):
    from calibration_process import operation as op

    cls = {"09": op.EC_operation_09, "12B": op.EC_operation_12B}[ver]
    return cls(
        tb_result_path=str(sp / tb_ref),
        fit_range=str(sp / "fit_range.json"),
        energy=str(sp / "ec_energy.json"),
        bkg_form=str(sp / "bkg_form.json"),
        x_path=str(DATA / ver / x_path),
        src_path=str(DATA / ver / src_path),
        save_path="x",
        save_fig_path="x",
        result_path="x",
    )


def legacy_file_config(ver, branch, mid):
    """Return the legacy ``file_config``-equivalent config tuple."""
    tb = legacy_tb(ver)
    ec = legacy_ec(ver)
    if branch == "tb":
        return tb.file_config(mid)
    if branch == "ec_source":
        return ec.src_config(mid)
    if branch == "ec_xray":
        return ec.xray_config(mid)
    raise ValueError(f"unknown branch {branch!r}")
