# -*- coding:utf-8 -*-
"""Timestamp helpers shared by the pipeline and plotting layers."""
import datetime
from pathlib import Path


def timestamp(format="%Y%m%d%H%M%S"):
    return datetime.datetime.now().strftime(format)


def headtime(path, forward=True):
    path = Path(path)
    if forward:
        new = path.parent / f"{timestamp()}_{path.stem}{path.suffix}"
    else:
        new = path.parent / f"{path.stem}_{timestamp()}{path.suffix}"
    return new
