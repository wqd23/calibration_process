# -*- coding:utf-8 -*-
"""
---------------------------------
INFO:

Module: GRID Plot Lib
Author: Wqd

GRID Collaboration/Calibration
---------------------------------
Features:

1. plot count-adc raw & fit graph
2. plot temp-bias fit graph
3. plot energy-channel fit graph
---------------------------------
"""

import grid_common.fit_utils as basic
from lib_reader.reader05.my_type import *
import nptyping as npt

from grid_common.naming import headtime
from grid_common.resolution import resolutionFunction

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os, math

# **************************************************************************************************************
# ******************************************* Raw/Fit plot function ********************************************
# **************************************************************************************************************


def raw_plot_single_channel(spectrum, x, title, ax=None, **kwargs):
    if ax is None:
        ax = plt.subplot()
    ax.step(x, spectrum, where="mid", label="raw data")
    ax.set_xlabel("ADC/channel")
    ax.set_ylabel("count rate/cps")
    ax.grid()
    ax.set_title(title)
    # kwargs
    # log scale
    if kwargs.get("log_scale", "linear") in ["log", "semilogx"]:
        ax.set_xscale("log")
    if kwargs.get("log_scale", "linear") in ["log", "semilogy"]:
        ax.set_yscale("log")
    return ax


def __raw_plot(
    spectrum: Float_array_4channel, x: Float_array_4channel, title: str, **kwargs
):
    default_kwargs = {"log_scale": "linear"}
    kwargs = {**default_kwargs, **kwargs}
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(4, 1, wspace=0.5, hspace=0.2, left=0.13, right=0.95)
    for ich, (spectrum_ich, x_ich) in enumerate(zip(spectrum, x)):
        ax = fig.add_subplot(gs[ich])
        if ich == 0:
            ax.set_title(title)
        # spectrum plot
        plt.step(x_ich, spectrum_ich, where="mid", label="raw data")
        ax.set_xlabel("ADC/channel")
        ax.set_ylabel("count rate/cps")
        ax.grid()
        ax.set_ylim([0, 1.2 * spectrum_ich.max()])
        # ax.set_xlim((0, 1000))
        # kwargs
        # log scale
        if kwargs["log_scale"] in ["log", "semilogx"]:
            ax.set_xscale("log")
        if kwargs["log_scale"] in ["log", "semilogy"]:
            ax.set_yscale("log")
        if kwargs.get("x_lim", None) is not None:
            ax.set_xlim(kwargs["x_lim"])
    return fig, fig.axes


def raw_plot(
    spectrum: Float_array_4channel,
    x: Float_array_4channel,
    title: str,
    save_path: Optional[str] = None,
    **kwargs,
) -> Optional[plt.Figure]:
    fig, axs = __raw_plot(spectrum, x, title, **kwargs)
    if save_path != None:
        fig.savefig(save_path)
        plt.close(fig)
    else:
        return fig


def fit_plot_single_channel(
    spectrum, x, fit_result, title, bkgForm, fit_range=None, ax=None, **kwargs
):
    if ax is None:
        ax = plt.subplot()
    if fit_result is None:
        ax.text(
            x[np.where(spectrum == np.max(spectrum))],
            max(spectrum),
            "None fit range offered\n no fit for this channel",
        )
        return ax
    if isinstance(fit_result, dict) and not fit_result.get("success", True):
        ax.text(
            x[np.where(spectrum == np.max(spectrum))],
            max(spectrum),
            f"fit FAILED:\n{fit_result.get('error', '')}",
        )
        return ax
    raw_plot_single_channel(spectrum, x, title, ax=ax, **kwargs)
    amplitude = fit_result["a"]
    center = fit_result["b"]
    sigma = fit_result["c"]
    fit_peak = basic.gaussianFunction([amplitude, center, sigma], x)
    bkg = fit_result["bkg"]
    fit_bkg = bkg["bkg_func"](x)
    fit_total = fit_peak + fit_bkg
    if fit_range is None:
        plot_range = [center - 3 * sigma, center + 3 * sigma]
    else:
        plot_range = fit_range
    qplot = (x >= plot_range[0]) * (x < plot_range[1])
    ax.plot(x[qplot], fit_peak[qplot], label="Gaussian Peak Fit")
    ax.plot(x[qplot], fit_bkg[qplot], label="Background Fit")
    ax.plot(x[qplot], fit_total[qplot], label="Gaussian Peak and Background Fit")
    info = f"""
    count rate = {fit_result["rate"]:.2e} $\\pm${fit_result["rate_err"]:.3e}
    center = {fit_result["b"]:.3e} $\\pm${fit_result["b_err"]:.3e}
    resolution = {fit_result["resolution"]:.3e} $\\pm${fit_result["resolution_err"]:.3e}
    """
    plt.text(
        center - 0.5 * sigma,
        max(fit_peak[qplot]) * 0.25,
        info,
        fontsize=10,
        bbox=dict(facecolor="pink", alpha=0.1),
        horizontalalignment="center",
        verticalalignment="center",
    )
    ax.legend(loc=0)
    if fit_range is not None:
        ax.axvline(plot_range[0])
        ax.axvline(plot_range[1])
    # set x,y lim
    ax.set_xlim(plot_range[0] - 3 * sigma, plot_range[1] + 3 * sigma)
    ax.set_ylim(0, 1.2 * np.max(spectrum[qplot]))
    return ax


def __fit_plot(
    spectrum: Float_array_4channel,
    x: Float_array_4channel,
    fit_result: List[Union[Dict[str, float], None]],
    title: str,
    bkgForm: str,
    fit_range=None,
    **kwargs,
):
    """universal plot function for peak fit

    Parameters
    ----------
    spectrum : Float_array_4channel
        spectrum count rate in x
    x : Float_array_4channel
        bin center
    fit_result : List[Union[Dict[str, float],None]]
        element is None or dict with keys [a,b,c,bkg]
    title : str
        fig title
    bkgForm : str
        'lin' or 'quad'
    fit_range :
        [[range_left, range_right]]*4
    save_path : str, optional
        full path to save the figure, by default None show the figure without save
    kwargs : Dict
        log_scale : str
            'linear', 'log', 'semilogx', 'semilogy'; by default 'linear'

    """
    default_kwargs = {"log_scale": "linear"}
    kwargs = {**default_kwargs, **kwargs}

    bkg_func = (
        lambda param, x: (param[0] * x + param[1])
        if bkgForm == "lin"
        else basic.quadFunction(param, x)
    )

    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(4, 1, wspace=0.5, hspace=0.2, left=0.13, right=0.95)
    for ich, (spectrum_ich, x_ich, fit_ich) in enumerate(zip(spectrum, x, fit_result)):
        ax = fig.add_subplot(gs[ich])
        if ich == 0:
            ax.set_title(title)
        # spectrum plot
        plt.step(x_ich, spectrum_ich, where="mid", label="raw data")
        if fit_ich != None:
            # fit plot data
            amplitude = fit_ich["a"]
            center = fit_ich["b"]
            sigma = fit_ich["c"]
            fit_peak = basic.gaussianFunction([amplitude, center, sigma], x_ich)
            bkg = fit_ich["bkg"]
            fit_bkg = bkg["bkg_func"](x_ich)
            fit_total = fit_peak + fit_bkg
            # fit plot
            if fit_range != None and fit_range[ich] != None:
                # fit range set
                plot_range = [fit_range[ich][0], fit_range[ich][1]]
            else:
                # fit range not set
                plot_range = [center - 3 * sigma, center + 3 * sigma]
            qplot = (x_ich >= plot_range[0]) * (x_ich < plot_range[1])
            plt.plot(x_ich[qplot], fit_peak[qplot], label="Gaussian Peak Fit")
            plt.plot(x_ich[qplot], fit_bkg[qplot], label="Background Fit")
            plt.plot(
                x_ich[qplot], fit_total[qplot], label="Gaussian Peak and Background Fit"
            )
            # fit info
            info = f"""
            count rate = {fit_ich["rate"]:.2e} $\\pm${fit_ich["rate_err"]:.3e}
            center = {fit_ich["b"]:.3e} $\\pm${fit_ich["b_err"]:.3e}
            resolution = {fit_ich["resolution"]:.3e} $\\pm${fit_ich["resolution_err"]:.3e}
            """
            plt.text(
                center - 0.5 * sigma,
                max(fit_peak[qplot]) * 0.25,
                info,
                fontsize=10,
                bbox=dict(facecolor="pink", alpha=0.1),
                horizontalalignment="center",
                verticalalignment="center",
            )
            # set x,y lim
            ax.set_xlim(plot_range[0] - 3 * sigma, plot_range[1] + 3 * sigma)
            ax.set_ylim(0, 1.2 * np.max(spectrum_ich[qplot]))
        else:
            plt.text(
                x_ich[np.where(spectrum_ich == np.max(spectrum_ich))][0],
                max(spectrum_ich),
                "None fit range offered\n no fit for this channel",
            )

        ax.set_xlabel("ADC/channel")
        ax.set_ylabel("count rante/cps")
        ax.legend(loc=0)
        ax.grid()
        # fit range
        if fit_range != None and fit_range[ich] != None:
            plt.axvline(fit_range[ich][0])
            plt.axvline(fit_range[ich][1])
        # kwargs
        # log scale
        if kwargs["log_scale"] in ["log", "semilogx"]:
            ax.set_xscale("log")
        if kwargs["log_scale"] in ["log", "semilogy"]:
            ax.set_yscale("log")
    return fig, fig.axes


def fit_plot(
    spectrum: Float_array_4channel,
    x: Float_array_4channel,
    fit_result: List[Union[Dict[str, float], None]],
    title: str,
    bkgForm: str,
    fit_range=None,
    save_path: Optional[str] = None,
    **kwargs,
):
    fig, axs = __fit_plot(spectrum, x, fit_result, title, bkgForm, fit_range, **kwargs)
    if save_path != None:
        fig.savefig(save_path)
        plt.close(fig)
    else:
        plt.show()


# ***************************************************************************************************************
# ******************************************* Temp-bias fit function ********************************************
# ***************************************************************************************************************


def plot_3d(xdata, ydata, zdata, xx, yy, zz, label, ax, title=None):
    ax.scatter3D(xdata, ydata, zdata, label="data point")  # type: ignore
    ax.plot_surface(xx, yy, zz, label="fit", alpha=0.4)  # type: ignore

    ax.set_title(title)
    ax.set_xlabel(label[0])
    ax.set_ylabel(label[1])
    ax.set_zlabel(label[2])  # type: ignore


def fit_err_plot_2d(
    data: npt.NDArray[npt.Shape["*,2"], npt.Floating],
    zdata: Float1D,
    func: Callable,
    label: Tuple[str, str, str],
    title: Optional[str] = None,
    save_path: Optional[str] = None,
):
    """plot for 2d fit

    Parameters
    ----------
    data : float array with shape [:,2]
        data[:,0] for xdata,data[:,1] for ydata
    zdata : Float1D
        data for z axis
    func : Callable
        func(x=data) return fit for zdata
    label : Tuple[str,str,str]
        x,y,z label string
    title : Optional[str], optional
        sup title of figure, by default None

    """
    xdata = data[:, 0]
    ydata = data[:, 1]
    x = np.arange(np.min(xdata) - 1, np.max(xdata) + 1)
    y = np.arange(np.min(ydata) - 1, np.max(ydata) + 1)
    xx, yy = np.meshgrid(x, y)
    zz = func(np.stack((xx, yy), axis=1))
    fig = plt.figure(figsize=(20, 10))
    if title != None:
        fig.suptitle(title)
    # 3d plot
    # fig_num = 2
    # ax = fig.add_subplot(fig_num,1,1,projection='3d')
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    plot_3d(xdata, ydata, zdata, xx, yy, zz, label, ax, title="3d plot")
    # 2d residual plot
    residual = func(data) - zdata
    ax = fig.add_subplot(1, 2, 2)
    cs = ax.tricontourf(xdata, ydata, residual / zdata * 100)  # type: ignore
    plt.colorbar(cs)
    ax.scatter(xdata, ydata, np.abs(residual / zdata * 100), c="r", marker="o")  # type: ignore
    ax.set_title("relative residual$\\times 100$")
    ax.set_xlabel(label[0])
    ax.set_ylabel(label[1])
    # set x,y lim to min-0.01(max-min), max+0.01(max-min), to make all data points in view
    ax.set_xlim(
        1.01 * np.min(xdata) - 0.01 * np.max(xdata),
        1.01 * np.max(xdata) - 0.01 * np.min(xdata),
    )
    ax.set_ylim(
        1.01 * np.min(ydata) - 0.01 * np.max(ydata),
        1.01 * np.max(ydata) - 0.01 * np.min(ydata),
    )
    plt.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)


# ********************************************************************************************************************
# ******************************************* Energy-channel fit function ********************************************
# ********************************************************************************************************************


def ec_plot(
    energy,
    center,
    result,
    src_energy,
    x_energy,
    src_result,
    x_result,
    save_path,
    energy_split_low,
    energy_split_high,
    per_channel_energy=None,
    xray_channels=None,
    src_energy_meas=None,
    x_energy_meas=None,
    ec_forms=None,
):
    x_energy, src_energy = np.asarray(x_energy), np.asarray(src_energy)
    # measurement-level energy arrays (one entry per measurement/group); they
    # can repeat (GRIDN1 measured Cs-137 twice and Co-60 on two files), unlike
    # the unique ``src_energy``/``x_energy`` above.
    src_e = np.asarray(src_energy if src_energy_meas is None else src_energy_meas)
    x_e = np.asarray(x_energy if x_energy_meas is None else x_energy_meas)
    # ``per_channel_energy`` lets channels carry different energy point sets
    # (GRIDN1: ch1/ch2 have X-ray + source, ch0/ch3 source only); ``xray_channels``
    # restricts the X-ray scatter to the channels that actually have those points.
    if per_channel_energy is None:
        en_by_ch = [np.asarray(energy)] * 4
    else:
        en_by_ch = [np.asarray(e) for e in per_channel_energy]
    center_low, center_high = [], []
    adc_low, adc_high = [], []
    energy_low, energy_high = [], []
    for _i in range(4):
        en_i = en_by_ch[_i]
        c = np.asarray(center[_i])
        cl = c[en_i < energy_split_low]
        chh = c[en_i > energy_split_high]
        center_low.append(cl)
        center_high.append(chh)
        lo = np.arange(np.min(cl), np.max(cl)) if cl.size else None
        hi = np.arange(np.min(chh), np.max(chh)) if chh.size else None
        adc_low.append(lo)
        adc_high.append(hi)
        energy_low.append(
            np.polyval(result[_i]["EC_low"], lo) if lo is not None else None
        )
        energy_high.append(
            np.polyval(result[_i]["EC_high"], hi) if hi is not None else None
        )

    def _col(groups, i, key):
        return np.array([
            g[i][key] if i < len(g) and g[i] is not None else np.nan for g in groups
        ])

    src_center = [_col(src_result, i, "b") for i in range(4)]
    src_center_err = [_col(src_result, i, "b_err") for i in range(4)]
    x_center = [_col(x_result, i, "b") for i in range(4)]
    x_center_err = [_col(x_result, i, "b_err") for i in range(4)]

    xpoint = (x_e > energy_split_high) | (x_e < energy_split_low)
    xpoint_not = np.logical_not(xpoint)
    gs = gridspec.GridSpec(
        2, 1, wspace=0.5, hspace=0.2, left=0.13, right=0.95, height_ratios=[4, 1]
    )
    for i in range(4):
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(gs[0])
        if xray_channels is None or i in xray_channels:
            ax.errorbar(
                x_center[i][xpoint],
                x_e[xpoint],
                xerr=x_center_err[i][xpoint],
                fmt="s",
                mfc="white",
                ms=6,
                elinewidth=1,
                capsize=3,
                barsabove=True,
                zorder=1,
                label=f" CH{i}",
            )
        ax.errorbar(
            src_center[i],
            src_e,
            xerr=src_center_err[i],
            fmt="^",
            mfc="white",
            ms=6,
            elinewidth=1,
            capsize=3,
            barsabove=True,
            zorder=0,
            label=f"source CH{i}",
        )
        if xray_channels is None or i in xray_channels:
            ax.errorbar(
                x_center[i][xpoint_not],
                x_e[xpoint_not],
                xerr=x_center_err[i][xpoint_not],
                fmt="s",
                mfc="red",
                ms=6,
                elinewidth=1,
                capsize=3,
                barsabove=True,
                zorder=1,
                label=f" CH{i} data not used",
            )

        if adc_low[i] is not None:
            form = "quadratic" if ec_forms is None else ec_forms[i]
            ax.plot(
                adc_low[i],
                energy_low[i],
                linestyle="-",
                label=f"{form} fit on EC data of ch{i}, < {energy_split_low}keV",
            )
        if adc_high[i] is not None:
            form = "quadratic" if ec_forms is None else ec_forms[i]
            ax.plot(
                adc_high[i],
                energy_high[i],
                linestyle="-",
                label=f"{form} fit on EC data of ch{i}, > {energy_split_high}keV",
            )
        ax.axhline(energy_split_low)
        ax.axhline(energy_split_high)
        ax.set_xlabel("ADC")
        ax.set_ylabel("Energy/keV")
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlim(80., 16384.0)
        ax.set_ylim(10.0, 1500.0)
        ax.legend(loc=0)
        ax.grid()
        fig.savefig(os.path.join(save_path, headtime(f"ec_fit_ch{i}.png")))

    e_union = np.concatenate((x_energy, src_energy))
    nonempty = [np.asarray(e) for e in en_by_ch if np.asarray(e).size]
    all_energy = np.concatenate(nonempty) if nonempty else np.asarray(energy)
    below = e_union[e_union < energy_split_low]
    above = e_union[e_union > energy_split_high]
    e_low = np.arange(np.min(all_energy), max(below)) if below.size else np.array([])
    e_high = np.arange(min(above), np.max(all_energy)) if above.size else np.array([])

    def _res_curve(e, coef):
        if coef is None or e.size == 0:
            return None
        return resolutionFunction(e, coef[0], coef[1], coef[2])

    resolution_low = [_res_curve(e_low, result[i]["resolution_low"]) for i in range(4)]
    resolution_high = [_res_curve(e_high, result[i]["resolution_high"]) for i in range(4)]
    src_resolution = [_col(src_result, i, "resolution") for i in range(4)]
    src_resolution_err = [_col(src_result, i, "resolution_err") for i in range(4)]
    x_resolution = [_col(x_result, i, "resolution") for i in range(4)]
    x_resolution_err = [_col(x_result, i, "resolution_err") for i in range(4)]

    gs = gridspec.GridSpec(
        2, 1, wspace=0.5, hspace=0.2, left=0.13, right=0.95, height_ratios=[4, 1]
    )
    for i in range(4):
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(gs[0])
        has_xray = xray_channels is None or i in xray_channels
        if has_xray and not any(map(math.isinf, x_resolution_err[i])):
            ax.errorbar(
                x_e[xpoint],
                x_resolution[i][xpoint] * 100,
                yerr=x_resolution_err[i][xpoint] * 100,
                fmt="s",
                mfc="white",
                ms=6,
                elinewidth=1,
                capsize=3,
                barsabove=True,
                zorder=1,
                label=f" CH{i} data used",
            )
            ax.errorbar(
                x_e[xpoint_not],
                x_resolution[i][xpoint_not] * 100,
                yerr=x_resolution_err[i][xpoint_not] * 100,
                fmt="s",
                mfc="white",
                ms=6,
                elinewidth=1,
                capsize=3,
                barsabove=True,
                zorder=1,
                label=f" CH{i} data not used",
            )
        elif has_xray:
            ax.scatter(x_e, x_resolution[i] * 100, label=f" CH{i}, inf in error")
        if src_resolution[i].size and not any(map(math.isinf, src_resolution_err[i])):
            ax.errorbar(
                src_e,
                src_resolution[i] * 100,
                yerr=src_resolution_err[i] * 100,
                fmt="^",
                mfc="white",
                ms=6,
                elinewidth=1,
                capsize=3,
                barsabove=True,
                zorder=0,
                label=f"source CH{i}",
            )
        elif src_resolution[i].size:
            ax.scatter(
                src_e, src_resolution[i] * 100, label=f"source CH{i}, inf in error"
            )
        if resolution_low[i] is not None:
            ax.plot(
                e_low,
                resolution_low[i] * 100,
                linestyle="-",
                label=f"Fit on resolution data of ch{i}, < {energy_split_low}keV",
            )
        if resolution_high[i] is not None:
            ax.plot(
                e_high,
                resolution_high[i] * 100,
                linestyle="-",
                label=f"Fit on resolution data of ch{i}, > {energy_split_high}keV",
            )
        ax.axvline(energy_split_low)
        ax.axvline(energy_split_high)
        ax.set_xlabel("Energy/keV")
        ax.set_ylabel("Resolution/%")
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlim(20, 1500)
        ymin = min(8, *list(src_resolution[i] * 100)) if src_resolution[i].size else 8
        ax.set_ylim(ymin, 150)
        ax.legend(loc=0)
        ax.grid()
        fig.savefig(os.path.join(save_path, headtime(f"resolution_fit_ch{i}.png")))
