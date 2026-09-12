# -*- coding:utf-8 -*-
"""Unit coverage for scripts/guess_fit_range.py (histogram peak guessing)."""

import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import guess_fit_range as g  # noqa: E402


def test_guess_histogram_synthetic_peak():
    rng = np.random.default_rng(0)
    center, sigma = 511.0, 30.0
    amp = rng.normal(center, sigma, 20000)
    res = g.guess_histogram(amp, bins=120, hist_range=(0, 1024))
    assert res["can_fit"], res["message"]
    lo, hi = res["fit_range"]
    assert lo < center < hi
    assert (hi - lo) > 4 * sigma          # sigma_multiple on each side


def test_guess_histogram_empty():
    res = g.guess_histogram([])
    assert res["can_fit"] is False
    assert res["fit_range"] is None


def test_guess_four_channels_and_cli(tmp_path):
    rng = np.random.default_rng(1)
    amps = {0: rng.normal(500, 25, 8000), 1: rng.normal(520, 28, 8000)}
    ranges = g.guess_four_channels(amps, channels=(0, 1))
    assert set(ranges) == {0, 1}

    npz = tmp_path / "a.npz"
    np.savez(npz, ch0=amps[0], ch1=amps[1])
    out = tmp_path / "fr.yaml"
    assert g.main([str(npz), "-o", str(out), "--measurement-id", "m"]) == 0
    d = yaml.safe_load(out.read_text())
    assert "m" in d and len(d["m"]) == 4
