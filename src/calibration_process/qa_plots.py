# -*- coding:utf-8 -*-
"""
Plot reduced-chi-square distributions per payload version and category (TB/EC).

Usage:
    python3 -m calibration_process.qa_plots [ver ...]

Saves PNG figures to data/qa_plots/<ver>_tb_hist.png and <ver>_ec_hist.png.
"""

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import util_lib as util
from .__init__ import CFG_PATH

cfg = util.load_config(CFG_PATH)

CHANNEL_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
CHANNEL_LABELS = ["ch0", "ch1", "ch2", "ch3"]


def _collect_redchi(root: Path):
    """Collect per-channel redchi values from fit pickles under root."""
    redchi = defaultdict(list)
    for pkl in sorted(root.glob("*.pickle")):
        try:
            data = util.pickle_load(pkl)
        except Exception:
            continue
        for ich, res in enumerate(data.get("fit_result", [])):
            if isinstance(res, dict) and res.get("success") is not False and "redchi" in res:
                redchi[ich].append(res["redchi"])
    return redchi


def _collect_tb2d_redchi(ver: str):
    """Collect per-channel redchi from TB 2D surface fit json."""
    best = None
    for f in sorted(Path(f"data/{ver}/tb_logs").glob("*_temp_bias_fit.json")):
        best = f
    if best is None:
        return {}
    channels = util.json_load(str(best))
    return {i: ch.get("redchi") for i, ch in enumerate(channels) if "redchi" in ch}


def _plot_hist(redchi_dict, title, outpath, tb2d_redchi=None):
    """Plot per-channel redchi histogram. If tb2d_redchi given, mark as vertical lines."""
    n_channels = len(redchi_dict)
    if n_channels == 0:
        return

    # Flatten to determine common bin range (exclude extreme outliers for binning)
    all_vals = []
    for vals in redchi_dict.values():
        all_vals.extend(vals)
    if not all_vals:
        return

    # Use log-spaced bins for wide-range distributions
    vmin = max(min(all_vals), 0.01)
    vmax = np.percentile(all_vals, 99.5) if len(all_vals) > 20 else max(all_vals)
    if vmax / vmin > 50:
        bins = np.geomspace(vmin, vmax * 1.2, 40)
        log_scale = True
    else:
        bins = np.linspace(vmin, vmax * 1.1, 40)
        log_scale = False

    fig, ax = plt.subplots(figsize=(8, 5))

    for ich in sorted(redchi_dict):
        vals = redchi_dict[ich]
        if not vals:
            continue
        ax.hist(
            vals,
            bins=bins,
            alpha=0.45,
            color=CHANNEL_COLORS[ich % 4],
            label=f"{CHANNEL_LABELS[ich]} (n={len(vals)}, p50={np.median(vals):.1f})",
            edgecolor="none",
        )

    if log_scale:
        ax.set_xscale("log")

    # Mark TB 2D redchi if provided
    if tb2d_redchi:
        for ich, val in tb2d_redchi.items():
            ax.axvline(
                val,
                color=CHANNEL_COLORS[ich % 4],
                linestyle="--",
                linewidth=1.2,
                alpha=0.7,
                label=f"{CHANNEL_LABELS[ich]} TB2D={val:.0f}" if ich < 4 else None,
            )

    ax.set_xlabel("reduced χ²")
    ax.set_ylabel("count")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  saved {outpath}")


def plot_version(ver: str):
    outdir = Path("data/qa_plots")
    outdir.mkdir(exist_ok=True)

    print(f"[{ver}]")

    # TB single-peak
    tb_redchi = _collect_redchi(Path(f"data/{ver}/single_process/TB_fit_result"))
    tb2d = _collect_tb2d_redchi(ver)
    _plot_hist(
        tb_redchi,
        f"{ver} — TB single-peak fit reduced χ²",
        outdir / f"{ver}_tb_hist.png",
        tb2d_redchi=tb2d,
    )

    # EC single-peak
    ec_redchi = _collect_redchi(Path(f"data/{ver}/single_process/EC_fit_result"))
    _plot_hist(
        ec_redchi,
        f"{ver} — EC single-peak fit reduced χ²",
        outdir / f"{ver}_ec_hist.png",
    )


def main():
    vers = sys.argv[1:] or list(cfg.keys())
    for ver in vers:
        plot_version(ver)
    print(f"\nAll plots saved to data/qa_plots/")


if __name__ == "__main__":
    main()
