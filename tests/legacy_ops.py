# -*- coding:utf-8 -*-
"""Factories that build the legacy operation objects for a version.

Used by the differential regression tests to compare the new explicit config
resolution against the historical ``file_config``-produced configs.  Paths and
the reader params are taken from the legacy ``config.json`` so these mirror the
historical constructor arguments exactly.
"""

from pathlib import Path

from calibration_process import operation as op
from calibration_process import util_lib as util
from calibration_process.__init__ import CFG_PATH

DATA = Path("data")


def _cfg(ver):
    return util.load_config(CFG_PATH)[ver]


def _sp(ver):
    return DATA / ver / "single_process"


_TBCLS = {
    "03B": op.TB_operation_03B, "04": op.TB_operation_04, "05B": op.TB_operation_05B,
    "07": op.TB_operation_07, "10B": op.TB_operation_10B, "11B": op.TB_operation_11B,
    "09": op.TB_operation_09, "12B": op.TB_operation_12B,
}
_ECCLS = {
    "03B": op.EC_operation_03B, "04": op.EC_operation_04, "05B": op.EC_operation_05B,
    "07": op.EC_operation_07, "10B": op.EC_operation_10B, "11B": op.EC_operation_11B,
    "09": op.EC_operation_09, "12B": op.EC_operation_12B,
}


def legacy_tb(ver):
    tb = _cfg(ver)["tb"]
    kwargs = dict(
        path=tb["path"],
        fit_range=str(_sp(ver) / "fit_range.json"),
        save_path=str(_sp(ver) / "TB_fit_result"),
        save_fig_path=str(_sp(ver) / "single_fit_fig"),
        result_path=str(DATA / ver / "tb_logs"),
    )
    if "file_map" in tb:
        kwargs["file_map"] = tb["file_map"]
    return _TBCLS[ver](**kwargs)


def legacy_ec(ver):
    ec = _cfg(ver)["ec"]
    base = dict(
        tb_result_path=ec["tb_result_path"],
        fit_range=str(_sp(ver) / "fit_range.json"),
        energy=str(_sp(ver) / "ec_energy.json"),
        bkg_form=str(_sp(ver) / "bkg_form.json"),
        x_path=ec["x_path"],
        src_path=ec["src_path"],
        save_path="x",
        save_fig_path="x",
        result_path="x",
    )
    if ver == "05B":  # base-class signature carries time_cut + x_config
        base.update(time_cut=ec["time_cut"], x_config=ec["x_config"])
    return _ECCLS[ver](**base)


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
