# -*- coding:utf-8 -*-
"""
process CLI interface for each GRID payload
----------
"""

import fire

from .cmd import VersionProcessOp, VersionProcessOp10B, VersionProcessOp12B
from . import util_lib as util
from .__init__ import CFG_PATH
from . import operation as op

cfg = util.load_config(CFG_PATH)

# (wrapper class, tb operation class, ec operation class, extra kwargs)
OPERATION_SPEC = {
    "03B": (VersionProcessOp, op.TB_operation_03B, op.EC_operation_03B,
            {"fp_method": "03"}),
    "04": (VersionProcessOp, op.TB_operation_04, op.EC_operation_04,
           {"fp_method": "04", "suffix": "txt"}),
    "05B": (VersionProcessOp, op.TB_operation_05B, op.EC_operation_05B,
            {"fp_method": None}),
    "07": (VersionProcessOp, op.TB_operation_07, op.EC_operation_07,
           {"fp_method": "07", "suffix": "txt"}),
    "10B": (VersionProcessOp10B, op.TB_operation_10B, op.EC_operation_10B,
            {"fp_method": "10"}),
    "11B": (VersionProcessOp10B, op.TB_operation_11B, op.EC_operation_11B,
            {"fp_method": "11"}),
    "09": (VersionProcessOp, op.TB_operation_09, op.EC_operation_09,
           {"fp_method": "09", "suffix": "txt"}),
    "12B": (VersionProcessOp12B, op.TB_operation_12B, op.EC_operation_12B,
            {"fp_method": "12"}),
}


def build_process(ver: str) -> VersionProcessOp:
    """Instantiate the process operation for a single payload version.

    Operations read data directories on construction, so they are only
    created here on demand instead of at import time.
    """
    wrapper_cls, tb_cls, ec_cls, kwargs = OPERATION_SPEC[ver]
    return wrapper_cls(
        tb_cls(**cfg[ver]["tb"]),
        ec_cls(**cfg[ver]["ec"]),
        **kwargs,
    )


class _LazyProcess:
    """Proxy that only builds the process operation on first attribute access."""

    def __init__(self, ver: str) -> None:
        self._ver = ver
        self._instance = None

    def _get_instance(self) -> VersionProcessOp:
        if self._instance is None:
            self._instance = build_process(self._ver)
        return self._instance

    def __getattr__(self, name: str):
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._get_instance(), name)

    def __dir__(self):
        # fire discovers members via dir(); trigger instantiation lazily here
        return sorted(set(super().__dir__()) | set(dir(self._get_instance())))


if __name__ == "__main__":
    fire.Fire({ver: _LazyProcess(ver) for ver in OPERATION_SPEC})
