# -*- coding:utf-8 -*-
"""
TB: temperature-bias heatmap of reduced chi-square (per channel)
EC: energy vs reduced chi-square scatter (per channel)

Usage:
    python3 -m calibration_process.qa_condition_plots [ver ...]
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

cfg = util.json_load(CFG_PATH)
CH_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def _tb_condition_plot(ver: str):
    """Heatmap: temp (x) vs bias (y), color = mean redchi across channels."""
    root = Path(f"data/{ver}/single_process/TB_fit_result")
    temps, biases, redchis = [], [], []
    for pkl in sorted(root.glob("*.pickle")):
        try:
            data = util.pickle_load(pkl)
        except Exception:
            continue
        fit = data.get("fit_result", [])
        tel = data.get("tel", {})
        # collect per-channel redchi that are valid
        vals = []
        for res in fit:
            if isinstance(res, dict) and res.get("success") is not False and "redchi" in res:
                vals.append(res["redchi"])
        if not vals:
            continue
        temp = np.mean(tel.get("tempSipm", [[np.nan]])[0])
        bias = np.mean(tel.get("bias", [[np.nan]])[0])
        temps.append(temp)
        biases.append(bias)
        redchis.append(np.mean(vals))

    if not temps:
        print(f"  [{ver}] TB: no valid data")
        return

    temps = np.array(temps)
    biases = np.array(biases)
    redchis = np.array(redchis)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharex=True, sharey=True)
    fig.suptitle(f"{ver} — TB reduced χ² vs temperature & bias", fontsize=13)

    for ich in range(4):
        ax = axes[ich]
        # collect per-channel redchi for this channel
        ch_temps, ch_biases, ch_redchis = [], [], []
        for pkl in sorted(root.glob("*.pickle")):
            try:
                data = util.pickle_load(pkl)
            except Exception:
                continue
            fit = data.get("fit_result", [])
            tel = data.get("tel", {})
            if ich >= len(fit):
                continue
            res = fit[ich]
            if not isinstance(res, dict) or res.get("success") is False or "redchi" not in res:
                continue
            temp = np.mean(tel.get("tempSipm", [[np.nan]])[0])
            bias = np.mean(tel.get("bias", [[np.nan]])[0])
            ch_temps.append(temp)
            ch_biases.append(bias)
            ch_redchis.append(res["redchi"])

        if not ch_temps:
            ax.set_title(f"ch{ich} (no data)")
            continue

        ch_temps = np.array(ch_temps)
        ch_biases = np.array(ch_biases)
        ch_redchis = np.array(ch_redchis)

        # cap redchi for color scale to avoid outliers dominating
        vmax = np.percentile(ch_redchis, 95) if len(ch_redchis) > 5 else ch_redchis.max()
        vmax = max(vmax, 1.5)  # at least show up to 1.5

        sc = ax.scatter(
            ch_temps, ch_biases,
            c=ch_redchis, cmap="RdYlGn_r",
            vmin=0.8, vmax=vmax,
            s=25, edgecolors="none",
        )
        ax.set_title(f"ch{ich} (n={len(ch_temps)})")
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Bias (V)")

    fig.colorbar(sc, ax=axes, label="reduced χ²", shrink=0.8)
    fig.tight_layout()
    outdir = Path("data/qa_plots")
    outdir.mkdir(exist_ok=True)
    outpath = outdir / f"{ver}_tb_heatmap.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  saved {outpath}")


def _ec_energy_plot(ver: str):
    """Scatter: energy (x) vs reduced chi-square (y), per channel."""
    root = Path(f"data/{ver}/single_process/EC_fit_result")
    energy_cfg = cfg[ver]["ec"]["energy"]
    energy_map = util.json_load(energy_cfg)

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(f"{ver} — EC reduced χ² vs photon energy", fontsize=13)

    for ich in range(4):
        energies, redchis = [], []
        for pkl in sorted(root.glob("*.pickle")):
            try:
                data = util.pickle_load(pkl)
            except Exception:
                continue
            fit = data.get("fit_result", [])
            fname = pkl.stem  # filename without .pickle
            if fname not in energy_map:
                continue
            e = energy_map[fname]
            if isinstance(e, list):
                e = e[ich] if ich < len(e) else np.nan
            if ich >= len(fit):
                continue
            res = fit[ich]
            if not isinstance(res, dict) or res.get("success") is False or "redchi" not in res:
                continue
            energies.append(e)
            redchis.append(res["redchi"])

        if not energies:
            continue

        energies = np.array(energies)
        redchis = np.array(redchis)

        ax.scatter(
            energies, redchis,
            color=CH_COLORS[ich], label=f"ch{ich}",
            s=30, alpha=0.7, edgecolors="none",
        )

    ax.set_xlabel("Photon energy (keV)")
    ax.set_ylabel("reduced χ²")
    ax.set_yscale("log")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    outdir = Path("data/qa_plots")
    outdir.mkdir(exist_ok=True)
    outpath = outdir / f"{ver}_ec_energy.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  saved {outpath}")


def main():
    vers = sys.argv[1:] or list(cfg.keys())
    for ver in vers:
        print(f"[{ver}]")
        _tb_condition_plot(ver)
        _ec_energy_plot(ver)
    print(f"\nAll plots saved to data/qa_plots/")


if __name__ == "__main__":
    main()
