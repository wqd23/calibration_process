# -*- coding:utf-8 -*-
"""Read-only QA analysis of persisted single-fit results.

For every ``*.fit.json`` under a version's ``single_process/*_fit_result``
directories, compute per-channel quality metrics (success, qa_flag, reduced
chi-square, boundary hits, peak position inside the fit window, relative center
error, resolution, residuals) and classify each channel as ``ok`` / ``suspect``
/ ``bad`` / ``fail`` / ``empty``.  Writes a machine-readable summary plus one
annotated PNG per measurement under ``data/{ver}/qa/`` (gitignored).

    python scripts/fit_qa.py 12B
    python scripts/fit_qa.py GRIDN1/EC --branch ec_xray
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from calibration_process import pipeline

BRANCH_DIR = {
    "tb": "TB_fit_result",
    "ec_source": "EC_fit_result",
    "ec_xray": "EC_fit_result",
    "neutron": "NEUTRON_fit_result",
}
DEFAULT_BRANCHES = ("tb", "ec_source", "ec_xray", "neutron")


def _load_spectrum(path: Path):
    if not path.exists():
        return None
    import polars as pl

    df = pl.read_parquet(path)
    out = {}
    for ch in df["channel"].unique().to_list():
        sub = df.filter(pl.col("channel") == ch).sort("bin")
        out[int(ch)] = (
            np.asarray(sub["x"], dtype=float),
            np.asarray(sub["spectrum"], dtype=float),
            np.asarray(sub["spectrum_err"], dtype=float),
        )
    return out


def _model(x, fit):
    a, b, c = fit["a"], fit["b"], fit["c"]
    if c <= 0:
        return np.zeros_like(x)
    peak = a * np.exp(-((x - b) ** 2) / (2 * c * c)) / (math.sqrt(2 * math.pi) * c)
    if fit.get("bkg", {}).get("bkg_info") == "lin" and "bk_slope" in fit:
        peak = peak + fit["bk_slope"] * x + fit["bk_intercept"]
    return peak


def _metrics(fit, spec, fit_range):
    m = {"success": bool(fit.get("success", False)), "qa_flag": fit.get("qa_flag")}
    m["redchi"] = float(fit.get("redchi", float("nan")))
    m["ndf"] = fit.get("ndf")
    m["boundary_hit"] = list(fit.get("boundary_hit", []))
    b, c = float(fit["b"]), float(fit["c"])
    b_err = float(fit.get("b_err", float("nan")))
    m["center"] = b
    m["sigma"] = c
    m["center_err"] = b_err
    m["resolution"] = float(fit.get("resolution", float("nan")))
    m["rel_center_err"] = abs(b_err / b) if b else float("nan")
    if fit_range and fit_range[0] is not None:
        lo, hi = float(fit_range[0]), float(fit_range[1])
        m["center_window_pos"] = (b - lo) / (hi - lo) if hi > lo else float("nan")
    else:
        m["center_window_pos"] = float("nan")
    if spec is not None:
        x, y, ye = spec
        if fit_range and fit_range[0] is not None:
            q = (x >= fit_range[0]) & (x <= fit_range[1])
            x, y, ye = x[q], y[q], ye[q]
        model = _model(x, fit)
        pull = (y - model) / np.where(ye > 0, ye, np.nan)
        m["residual_rms"] = float(np.sqrt(np.nanmean((y - model) ** 2)))
        m["max_pull"] = float(np.nanmax(np.abs(pull))) if pull.size else float("nan")
    else:
        m["residual_rms"] = float("nan")
        m["max_pull"] = float("nan")
    m["category"] = _classify(m)
    return m


def _classify(m):
    if not m["success"] or m["qa_flag"] == "fail":
        return "fail"
    pos = m.get("center_window_pos", float("nan"))
    if m["boundary_hit"]:
        return "bad"
    if not math.isnan(pos) and (pos < 0.08 or pos > 0.92):
        return "bad"
    if not math.isnan(m["rel_center_err"]) and m["rel_center_err"] > 0.05:
        return "bad"
    if not math.isnan(m["max_pull"]) and m["max_pull"] > 8:
        return "suspect"
    if not math.isnan(m["redchi"]) and m["redchi"] > 5:
        return "suspect"
    return "ok"


def _plot(out, branch, mid, fits, spec, ranges):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 4, figsize=(18, 7), sharex=False)
    for ch in range(4):
        fit = fits[ch] if ch < len(fits) else None
        ax, rax = axes[0][ch], axes[1][ch]
        ax.set_title(f"{mid} ch{ch}")
        if spec is not None and ch in spec:
            x, y, ye = spec[ch]
            ax.errorbar(x, y, yerr=ye, fmt=".", ms=2, lw=0.5, color="tab:blue")
            rng = ranges[ch] if ch < len(ranges) else None
            if rng and rng[0] is not None:
                ax.axvline(rng[0], color="k", lw=0.5, ls="--")
                ax.axvline(rng[1], color="k", lw=0.5, ls="--")
            if fit:
                q = np.ones_like(x, dtype=bool)
                if rng and rng[0] is not None:
                    q = (x >= rng[0]) & (x <= rng[1])
                xx = x[q]
                ax.plot(xx, _model(xx, fit), "r-", lw=1.0)
                rax.plot(xx, (y[q] - _model(xx, fit)) / np.where(ye[q] > 0, ye[q], np.nan), ".-", ms=2)
        rax.axhline(0, color="k", lw=0.4)
        if fit:
            m = _metrics(fit, spec[ch] if spec and ch in spec else None,
                         ranges[ch] if ch < len(ranges) else None)
            rax.set_title(
                f"{m['category']} redchi={m['redchi']:.1f} c={m['center']:.1f}", fontsize=8
            )
        else:
            rax.set_title("empty", fontsize=8)
    fig.suptitle(branch)
    fig.tight_layout()
    fig.savefig(out / f"{branch}_{mid}.png", dpi=80)
    plt.close(fig)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("version")
    p.add_argument("--branch", action="append", choices=DEFAULT_BRANCHES)
    p.add_argument("--out", default=None)
    args = p.parse_args(argv)

    ver = args.version
    rt = pipeline.load_rt(ver)
    out = Path(args.out) if args.out else rt.data_dir / "qa"
    fig_dir = out / "fig"
    fig_dir.mkdir(parents=True, exist_ok=True)
    branches = args.branch or [b for b in DEFAULT_BRANCHES
                               if (rt.data_dir / "single_process" / BRANCH_DIR[b]).is_dir()]

    summary = {}
    rows = []
    for branch in branches:
        sub = rt.data_dir / "single_process" / BRANCH_DIR[branch]
        if not sub.is_dir():
            continue
        for fj in sorted(sub.glob("*.fit.json")):
            mid = fj.name[: -len(".fit.json")]
            data = json.loads(fj.read_text())
            fits = data["fit_result"]
            try:
                ranges = rt.fit_range(branch, mid)
            except KeyError:
                continue
            spec = _load_spectrum(sub / f"{mid}.spectrum.parquet")
            entry = {}
            for ch in range(4):
                if ch >= len(fits) or fits[ch] is None:
                    entry[str(ch)] = {"category": "empty"}
                    continue
                cspec = spec[ch] if spec and ch in spec else None
                entry[str(ch)] = _metrics(fits[ch], cspec, ranges[ch] if ch < len(ranges) else None)
                m = entry[str(ch)]
                rows.append([ver, branch, mid, ch, m["category"], m["qa_flag"],
                             m["center"], m["sigma"], m["redchi"], m["rel_center_err"],
                             m["center_window_pos"], ";".join(m["boundary_hit"])])
            summary.setdefault(branch, {})[mid] = entry
            _plot(fig_dir, branch, mid, fits, spec, ranges)

    (out / "qa_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    with (out / "qa_summary.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ver", "branch", "id", "ch", "category", "qa_flag", "center",
                    "sigma", "redchi", "rel_center_err", "center_window_pos", "boundary_hit"])
        w.writerows(rows)
    print(f"wrote {out}/qa_summary.json and {len(rows)} channel rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
