# -*- coding:utf-8 -*-
"""Row-wise sugar for the columnar L1 frames passed to selection callbacks."""
from __future__ import annotations

import numpy as np


class EventTable:
    """Thin per-row accessor over a ``{column: 1-D-or-2-D array}`` frame table.

    Selection callbacks receive the raw frame mapping directly (columnar, so
    whole-array cuts stay vectorized); ``EventTable`` is only for the occasional
    per-event inspection, e.g. ``EventTable(frames).row(0)["waveform_data"]``.
    """

    def __init__(self, frames):
        self._frames = frames

    def __len__(self):
        return len(next(iter(self._frames.values())))

    def column(self, name):
        return self._frames[name]

    def row(self, i):
        return {k: np.asarray(v)[i] for k, v in self._frames.items()}
