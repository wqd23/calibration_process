# -*- coding:utf-8 -*-
"""Byte-source helpers for the payload readers.

* binary files (10B/11B/12B, 03B/05B) are read as a ``uint8`` stream with
  :func:`load_binary`;
* space-separated ASCII hex text files (04/07/09) are decoded with
  :func:`load_hex_text`.

The returned byte stream is consumed by :mod:`packet_parser` (A/C) or by the
03B/05B decoders in :mod:`reader05.readout`.
"""
import numpy as np


def load_binary(path) -> np.ndarray:
    """Read a binary file as a ``uint8`` byte stream."""
    return np.fromfile(path, dtype=np.uint8)


def load_hex_text(path) -> np.ndarray:
    """Decode a space-separated hex-text file (04/07/09) to a byte stream."""
    raw = np.fromfile(path, dtype=np.uint8).tobytes()
    text = raw.decode("ascii", errors="ignore")
    vals = []
    for tok in text.split():
        try:
            vals.append(int(tok, 16))
        except ValueError:
            pass
    return np.asarray(vals, dtype=np.uint8)
