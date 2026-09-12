# -*- coding:utf-8 -*-
"""Suggest per-channel Gaussian fit ranges from an amplitude histogram.

Ported from ``gridN_cali/gaus_guess.py`` (the GRID-N1 reference): histogram the
amplitudes, let lmfit guess a Gaussian, run a trial fit, gate on R^2 and the
relative sigma error, then set ``fit_range = center +/- sigma_multiple * sigma``.
The output is directly usable as a ``fit_range_*.yaml`` channel entry.

``amp`` values are in ADC (the reader's ``data_max - data_base/4``).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml
from lmfit.models import GaussianModel


def _fit(counts, centers, bin_width, model, min_r2, max_sigma_stderr_ratio, sigma_multiple):
    lm = GaussianModel()
    try:
        params = lm.guess(counts, x=centers)
    except Exception as exc:  # noqa: BLE001
        return {"can_fit": False, "message": f"guess failed: {exc}", "fit_range": None, "p0": None}

    x_lo, x_hi = float(centers.min()), float(centers.max())
    x_span = max(x_hi - x_lo, 1.0)
    params["center"].set(min=max(x_lo - 0.3 * x_span, 0.0), max=x_hi + 0.3 * x_span)
    params["sigma"].set(min=max(bin_width * 0.5, 1.0), max=x_span * 3)
    try:
        trial = lm.fit(counts, params, x=centers)
    except Exception as exc:  # noqa: BLE001
        return {"can_fit": False, "message": f"trial fit failed: {exc}", "fit_range": None, "p0": None}

    fitted = trial.best_fit
    ss_res = float(np.sum((counts - fitted) ** 2)) if fitted is not None else np.nan
    ss_tot = float(np.sum((counts - np.mean(counts)) ** 2)) if counts.size > 1 else np.nan
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    center = float(trial.params["center"].value)
    sigma = abs(float(trial.params["sigma"].value))
    err = trial.params["sigma"].stderr

    fail = None
    if not np.isfinite(r2) or r2 < min_r2:
        fail = f"r2={r2:.3f} < {min_r2}"
    elif center < 0.0 or center > x_hi:
        fail = f"center {center:.1f} outside [0, {x_hi:.0f}]"
    elif err is not None and sigma > 0 and err / sigma > max_sigma_stderr_ratio:
        fail = f"sigma relative error {err / sigma:.3f} > {max_sigma_stderr_ratio}"
    if fail is not None:
        return {"can_fit": False, "message": fail, "fit_range": None, "p0": None}

    fit_range = [max(0.0, center - sigma_multiple * sigma),
                 center + sigma_multiple * sigma]
    return {
        "can_fit": True,
        "message": "ok",
        "fit_range": [float(fit_range[0]), float(fit_range[1])],
        "p0": {"amplitude": float(trial.params["amplitude"].value),
               "center": center, "sigma": sigma},
    }


def guess_histogram(amp, *, bins=100, hist_range=None, min_r2=0.85,
                    max_sigma_stderr_ratio=0.5, sigma_multiple=4.0):
    """Guess ``[lo, hi]`` for one channel, or a dict with ``can_fit=False``."""
    amp = np.asarray(amp, dtype=float).reshape(-1)
    if amp.size == 0:
        return {"can_fit": False, "message": "no data", "fit_range": None, "p0": None}
    if hist_range is None:
        lo, hi = float(np.min(amp)), float(np.max(amp))
        if hi <= lo:
            hi = lo + 1.0
        hist_range = (lo, hi)
    counts, edges = np.histogram(amp, bins=bins, range=hist_range)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = float(np.diff(centers).mean()) if centers.size > 1 else 1.0
    return _fit(counts, centers, bin_width, "gauss", float(min_r2),
                float(max_sigma_stderr_ratio), float(sigma_multiple))


def guess_four_channels(amps, *, channels=(0, 1, 2, 3), **kwargs):
    """Guess a fit range for each channel; returns ``{channel: range|None}``."""
    return {int(ch): guess_histogram(amps[int(ch)], **kwargs)["fit_range"]
            for ch in channels}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("npz", help="npz with one amplitude array per channel (ch0..ch3)")
    p.add_argument("-o", "--out", help="write {channel: range} as YAML")
    p.add_argument("--measurement-id", help="wrap the YAML as a fit_range measurements entry")
    p.add_argument("--bins", type=int, default=100)
    args = p.parse_args(argv)

    data = np.load(args.npz)
    amps = {int(k[2:]): data[k] for k in data.files if k.startswith("ch")}
    ranges = guess_four_channels(amps, channels=tuple(sorted(amps)))
    chans = [ranges.get(ch) for ch in range(4)]
    print(yaml.safe_dump({"channels": ranges}, sort_keys=True).strip())
    if args.out:
        payload = ({args.measurement_id: chans} if args.measurement_id
                   else {"channels": ranges})
        Path(args.out).write_text(yaml.safe_dump(payload, sort_keys=False))
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
