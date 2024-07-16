# -*- coding:utf-8 -*-
"""
process CLI interface for each GRID payload
----------
"""
import fire

from .cmd import VersionProcessOp
from . import util_lib as util
from .__init__ import CFG_PATH
from . import operation as op

cfg = util.json_load(CFG_PATH)

process03B = VersionProcessOp(
    op.TB_operation_03B(**cfg["03B"]["tb"]),
    op.EC_operation_03B(**cfg["03B"]["ec"]),
    fp_method="03",
)
process04 = VersionProcessOp(
    op.TB_operation_04(**cfg["04"]["tb"]),
    op.EC_operation_04(**cfg["04"]["ec"]),
    fp_method="04",
)
process05B = VersionProcessOp(
    op.TB_operation_05B(**cfg["05B"]["tb"]),
    op.EC_operation_05B(**cfg["05B"]["ec"]),
    fp_method=None,
)

process07 = VersionProcessOp(
    op.TB_operation_07(**cfg["07"]["tb"]),
    op.EC_operation_07(**cfg["07"]["ec"]),
    fp_method="07",
)

process10B = VersionProcessOp(
    op.TB_operation_10B(**cfg["10B"]["tb"]),
    op.EC_operation_10B(**cfg["10B"]["ec"]),
    fp_method="10",
)

if __name__ == "__main__":
    fire.Fire({
        "03B": process03B,
        "04": process04,
        "05B": process05B,
        "07": process07,
        "10B": process10B
    })
