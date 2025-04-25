# -*- coding:utf-8 -*-
"""
GRID version independent Functions Library
----------

Basic functions for GRID data processing\n
"""
import lib_reader.reader05.gridBasicFunctions as basic
from lib_reader.reader05.my_type import *

import numpy as np
import lmfit.models as fitmodel
import lmfit
from scipy.optimize import curve_fit
import dill as pickle
import json
import os
import datetime
import re
from pathlib import Path


def count_spectrum(amp, nbins, spec_range, bin_width, adc_max):
    spectrum, x = basic.getSpectrum(
        amp, nbins=nbins, specRange=spec_range, binWidth=bin_width, adcMax=adc_max
    )
    spectrum_err: Float_array_4channel = [
        basic.gehrelsErr(spectrum_ich) for spectrum_ich in spectrum
    ]  # type: ignore
    return spectrum, spectrum_err, x


def count(
    amp: Float_array_4channel, bin_width: float, spec_range: List[List[float]]
) -> Tuple[Float_array_4channel, Float_array_4channel, Float_array_4channel]:
    spectrum_all = []
    spectrum_err_all = []
    x_all = []
    for a, s_range in zip(amp, spec_range):
        binEdges = np.arange(s_range[0], s_range[1] + bin_width, bin_width)
        spectrum, x = np.histogram(a, bins=binEdges)
        # assert spectrum.sum() == a.size, "spectrum sum is not equal to the number of events"
        x = (x[:-1] + x[1:]) / 2
        spectrum_err = basic.gehrelsErr(spectrum)
        spectrum_all.append(spectrum)
        spectrum_err_all.append(spectrum_err)
        x_all.append(x)
    return spectrum_all, spectrum_err_all, x_all


class FitError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
        # message for print
        self.message = args[-1]


def peak_fit(
    x: Float1D, data: Float1D, error: Float1D, bkg_form: Optional[str] = None
) -> Dict[str, Any]:
    """peak fit for single channel

    Parameters
    ----------
    x : Float1D
        bin center
    data : Float1D
        count rate
    error : Float1D
        data error
    bkg_form : Optional[str], optional
        'lin','quad', 'exp', 'guas' by default None

    Returns
    -------
    Dict[str, Any]
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'peak_amplitude': area of the gaussian peak,
            'peak_center': center of the gaussian peak,
            'peak_sigma': standard deviation of the gaussian peak,
            'peak_amplitude_err': error of area of the gaussian peak,
            'peak_center_err': error of center of the gaussian peak,
            'peak_sigma_err': error of standard deviation of the gaussian peak,
            'bkg': {
                'bkg_info': 'lin','quad', 'exp' or None
                'bkg_func': function with bkg best value
                **bkg_parameter
            }
        }\n
    with fit function as\n
    `y =  peak_amplitude * exp(-(x - peak_center) ** 2 / (2 * peak_sigma ** 2)) / sqrt(2 * pi * peak_sigma ^ 2) + bk_a * x  + bk_b`\n

    Raises
    -------
    FitError : when fit is not success or uncertainties could not be estimated
       result : lmfit.model.ModelResult
       des : str
            description of the FitError
    """
    mod = fitmodel.GaussianModel(prefix="peak_")
    param = mod.guess(data, x=x)
    if bkg_form == "lin":
        b_mod = lmfit.models.LinearModel(prefix="bk_")
    elif bkg_form == "quad":
        b_mod = lmfit.models.QuadraticModel(prefix="bk_")
    elif bkg_form == "exp":
        b_mod = lmfit.models.ExponentialModel(prefix="bk_")
    else:
        b_mod = fitmodel.GaussianModel(prefix="bk_")
    if bkg_form != None:
        # while use bkg, right half to guess data
        data_param = mod.guess(data[len(data) // 2 :], x=x[len(data) // 2 :])
        bkg_data = data - mod.eval(data_param, x=x)
        b_param = b_mod.guess(bkg_data, x=x)
        mod = mod + b_mod
        param = data_param + b_param
    # set bound for param
    # assume center is in the xdata range, xdata range bigger than 3*sigma(half the peak)
    # param['peak_amplitude'].min, param['peak_amplitude'].max = 0, 10*np.max(ydata)
    param["peak_amplitude"].min = 0
    param["peak_center"].min, param["peak_center"].max = x[0], x[-1]
    param["peak_sigma"].min, param["peak_sigma"].max = 0, (x[-1] - x[0]) / 3.0
    param.add("peak_resolution", expr="peak_fwhm/peak_center")
    result = mod.fit(data, param, x=x, weights=1.0 / error)
    fitResult = {k: v for k, v in result.values.items() if "bk_" not in k}
    fitResult.update(
        {
            "peak_amplitude_err": result.params["peak_amplitude"].stderr,
            "peak_center_err": result.params["peak_center"].stderr,
            "peak_sigma_err": result.params["peak_sigma"].stderr,
        }
    )

    if bkg_form == "lin":
        fitResult.update(
            {
                "bk_a": result.params["bk_slope"].value,
                "bk_b": result.params["bk_intercept"].value,
            }
        )
    if bkg_form != None:
        fitResult["bkg"] = {
            "bkg_info": bkg_form,
            "bkg_func": lambda x: b_mod.eval(x=x, **result.best_values),
        }
        fitResult["bkg"].update({k: v for k, v in result.values.items() if "bk_" in k})
    else:
        fitResult["bkg"] = {
            "bkg_info": "None",
            "bkg_func": lambda x: 0 * x,
        }
    # fit failed or can't get uncertainty
    if not result.success or any(map(lambda x: x == None, fitResult.values())):
        # param too much, that can't get uncertainty
        if bkg_form == "lin":
            # fit without bkg
            try:
                fitResult = peak_fit(x, data, error)
            except FitError as e:
                raise FitError(
                    e.args[0],
                    f"failed to fit with {bkg_form} background and failed to fit without background",
                )
            # add the linear zero bkg
            fitResult.update(
                {
                    "bk_a": 0,
                    "bk_b": 0,
                }
            )
        else:
            raise FitError(result, f"failed to fit with {bkg_form} background")
    return fitResult


def tempbias2DFunctionInternal(x, G0, k, V0, b, c):
    """
    Auxiliary function to calculate the value of 2D temperature-bias response function for correlative temperature-bias fit, used for ODR fit\n
    in the form of\n
    `f(x, y) = G0 * (y - k * x - V0) ^ 2 * (-x ^ 2 + b * x + c)`
    Parameters
    ----------
    param : list,
        parameters of the function, in the form of `[G0, k, V0, b, c]`\n
    x : array-like,
        array containing both x and y data, being `xdata = x[:,0]` and `ydata = x[:,1]`\n
    Returns
    ----------
    float or array-like,
        values of temperature-bias response function corresponding to input x and y\n
    """
    xdata = x[:, 0]
    ydata = x[:, 1]
    Vov = ydata - k * xdata - V0
    return G0 * Vov**2 * (-(xdata**2) + b * xdata + c)


def tempbias2DFunction(temp, bias, G0, k, V0, b, c):
    Vov = bias - k * temp - V0
    return G0 * Vov**2 * (-(temp**2) + b * temp + c)


def resolutionFunction(x, a, b, c):
    y = np.sqrt(a * x * x + b * x + c) / x
    y[a * x * x + b * x + c < 0] = 0.0
    return y


def temp_bias_lmfit(
    center: Float1D,
    center_err: Float1D,
    temp: Float1D,
    temp_err: Float1D,
    bias: Float1D,
    bias_err: Float1D,
):

    data = np.stack([temp, bias, center], axis=1)
    error = np.stack([temp_err, bias_err, center_err], axis=1)
    model = lmfit.Model(tempbias2DFunctionInternal)
    # param = model.guess(data[:,2], data[:,0:2])

    model.set_param_hint("G0", value=0.030317363114139174, min=1e-5, max=100.0)
    model.set_param_hint("k", value=18.9e-3, min=1e-3, max=50e-3)
    model.set_param_hint("V0", value=24.6, min=23.0, max=25.0)
    model.set_param_hint("b", value=0.0, min=-500.0, max=500.0)
    model.set_param_hint("c", value=1e4, min=1e3, max=1e8)
    param = model.make_params()
    result = model.fit(
        data[:, 2],
        x=data[:, 0:2],
        params=param,
        weights=1.0 / error[:, 2],
        method="trf",
    )
    fitResult = {
        "G0": result.best_values["G0"],
        "k": result.best_values["k"],
        "V0": result.best_values["V0"],
        "b": result.best_values["b"],
        "c": result.best_values["c"],
        "G0_err": param["G0"].stderr,
        "k_err": param["k"].stderr,
        "V0_err": param["V0"].stderr,
        "b_err": param["b"].stderr,
        "c_err": param["c"].stderr,
        "chisquare": result.chisqr,
    }
    if not result.success or any(map(lambda x: x == None, fitResult.values())):
        raise FitError(result, "faled to do temp bias 2d fit")


def temp_bias_fit_curvefit(
    center: Float1D, center_err: Float1D, temp: Float1D, bias: Float1D
) -> Dict[str, float]:
    data = np.stack([temp, bias], axis=1)
    initial_guess = [0.00078, 18.9e-3, 24.6, 0, 1e4]
    try:
        popt, pcov = curve_fit(
            tempbias2DFunctionInternal,
            data,
            center,
            sigma=center_err,
            absolute_sigma=True,
            p0=initial_guess,
        )
    except RuntimeError as e:
        raise FitError(e.args[0])
    perr = np.sqrt(np.diag(pcov))
    param = popt.tolist()
    chisq = sum(
        (basic.residualTempbias2D(param, temp, bias, center) / center_err) ** 2
    ) / (len(center_err) - len(param))
    fitResult = {
        "G0": param[0],
        "k": param[1],
        "V0": param[2],
        "b": param[3],
        "c": param[4],
        "G0_err": perr[0],
        "k_err": perr[1],
        "V0_err": perr[2],
        "b_err": perr[3],
        "c_err": perr[4],
        "chisquare": chisq,
    }
    return fitResult


def pickle_save(data, path: str):
    """save basic data with pickle

    Parameters
    ----------
    data : Any
        basic python data
    path : str
        full path to save the data
    """
    with open(path, "wb") as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def pickle_load(path) -> Any:
    """load data saved by pickle

    Parameters
    ----------
    path : str
        full path of data

    Returns
    -------
    Any
        data
    """
    with open(path, "rb") as f:
        return pickle.load(f)


def json_save(data, path: str):
    """save the data into a json file

    Parameters
    ----------
    data : Any
        any oject could save as json
    path : str
        full path "path/to/file.json"
    """
    with open(path, "w") as f:
        f.write(json.dumps(data))


def timestamp(format="%Y%m%d%H%M%S"):
    return datetime.datetime.now().strftime(format)


def json_time_save(data, path: str, forward=False):
    """save data with a time appended file name, wrapper for json_save

    Parameters
    ----------
    data : Any
        data object
    path : str
        full path, time info will be appended to it
    forward : bool
        add time stamp before
    """
    new = headtime(path, forward)
    json_save(data, new)


def json_headtime_save(data, path: str):
    json_time_save(data, path, forward=True)


def headtime(path, forward=True):
    path = Path(path)
    if forward:
        new = path.parent / f"{timestamp()}_{path.stem}{path.suffix}"
    else:
        new = path.parent / f"{path.stem}_{timestamp()}{path.suffix}"
    return new


def json_load(path: str):
    with open(path, "r") as f:
        return json.loads(f.read())


def tb_corr(
    corr: List[Callable[[float, float], float]],
    temp: Float_array_4channel,
    bias: Float_array_4channel,
) -> Float_4channel:
    temp_avg = [np.average(t) for t in temp]
    bias_avg = [np.average(b) for b in bias]
    corr_factor = [f(t, b) for f, t, b in zip(corr, temp_avg, bias_avg)]
    return corr_factor


def resolution_polyfit(energy: Float1D, resolution: Float1D, resolution_err: Float1D):
    data = (resolution * energy) ** 2
    error = 2 * energy**2 * resolution_err
    p0, pcov = np.polyfit(energy, data, deg=2, w=1.0 / error, cov=True)
    perr = np.sqrt(np.diag(pcov))
    popt = list(p0)
    perr = list(perr)
    return p0, pcov


def resolution_ExprFit(energy: Float1D, resolution: Float1D, resolution_err: Float1D):
    p0, pcov = resolution_polyfit(energy, resolution, resolution_err)
    mod = lmfit.models.ExpressionModel(
        expr="sqrt(a*x*x+b*x+c)/x", independent_vars=["x"]
    )
    param = mod.make_params()
    param["a"].min, param["b"].min, param["c"].min = 0, 0, 0
    xOr01 = lambda x: x if x > 0 else 0.1
    param["a"].value, param["b"].value, param["c"].value = (
        xOr01(p0[0]),
        xOr01(p0[1]),
        xOr01(p0[2]),
    )
    result = mod.fit(resolution, param, weights=1.0 / resolution_err, x=energy)
    popt = [result.best_values["a"], result.best_values["b"], result.best_values["c"]]
    perr = [param["a"].stderr, param["b"].stderr, param["c"].stderr]
    return popt, perr


def resolution_lmfit(energy: Float1D, resolution: Float1D, resolution_err: Float1D):
    # initial value
    p0, pcov = resolution_polyfit(energy, resolution, resolution_err)
    mod = lmfit.Model(resolutionFunction)
    param = mod.make_params(
        a=0.1 if p0[0] > 0 else p0[0],
        b=0.1 if p0[0] > 0 and p0[2] > 0 else p0[1],
        c=0.1 if p0[2] > 0 else p0[2],
    )
    # if p0[0]>0 and p0[2]>0:
    # param['a'].min=0
    param["b"].min = 0
    param["c"].min = 0
    result = mod.fit(resolution, param, weights=1.0 / resolution_err, x=energy)
    popt = [result.best_values["a"], result.best_values["b"], result.best_values["c"]]
    perr = [param["a"].stderr, param["b"].stderr, param["c"].stderr]
    return popt, perr


def get_key(file_name: str):
    assert os.path.exists(file_name), f"{file_name} does not exist!"
    return os.path.basename(file_name)


def get_hpge_data(path):
    data = np.loadtxt(path)
    counts = data
    # unit: keV
    energy = np.arange(0, len(counts), 1)
    return energy, counts


def get_fit_dict(path: str, ec_energy, suffix="dat"):
    """get energy and fit result of both x and src files, return result dicts

    Parameters
    ----------
    path : str
        path of EC_fit_result
    ec_energy : Any
        ec_energy through load_json
    suffix : str
        [Default: 'dat'] suffix of data files, 'txt' for 07 and 'dat' for 03 & 05
    """
    # key -> energy, value -> fit-result
    x_result = {}
    src_result = {}
    f_list = os.listdir(path)
    # X: for 03 & 07, [?p?.pickle] alike gauge
    regex_x_energy = r"^(\d+)p(\d+)(.pickle)$"
    # X: for 05, [X/x M/m···observe.pickle] alike gauge
    regex_xm = r"[xX][mM]\S+(observe.pickle)$"
    # SRC: [src···.pickle] alike gauge
    regex_src = r"src\S+(.pickle)$"
    for f in f_list:
        match_x_energy = re.search(regex_x_energy, f)
        match_xm = re.search(regex_xm, f)
        match_src = re.search(regex_src, f)
        if match_x_energy != None:
            integer, fraction, _ = match_x_energy.groups()
            # ec_energy x key: [?p?] alike gauge
            f_name = integer + "p" + fraction
            energy = ec_energy[f_name]
            data = pickle_load(os.path.join(path, f))
            x_result[energy] = data["fit_result"]
            continue
        elif match_xm != None:
            # ec_energy x key: [filename] alike gauge
            f_name = re.sub("pickle", suffix, f)
            energy = ec_energy[f_name]
            data = pickle_load(os.path.join(path, f))
            x_result[energy] = data["fit_result"]
            continue
        elif match_src != None:
            # ec_energy src key: [filename] alike gauge
            f_name = re.sub("pickle", suffix, f)
            energy = ec_energy[f_name]
            data = pickle_load(os.path.join(path, f))
            src_result[energy] = data["fit_result"]
    return x_result, src_result


def print_temp(tel):
    temp = tel["tempSipm"]
    avg = [np.mean(temp[i]) for i in range(4)]
    std = [np.std(temp[i]) for i in range(4)]
    print(
        f"CH0: {avg[0]:.2f}+/-{std[0]:.2f} CH1: {avg[1]:.2f}+/-{std[1]:.2f} CH2: {avg[2]:.2f}+/-{std[2]:.2f} CH3: {avg[3]:.2f}+/-{std[3]:.2f}"
    )


def data_save(sci, tel, path, name):
    dir = os.path.join(path, name)
    if not os.path.exists(dir):
        os.makedirs(dir)
    keys = [
        "amp",
        "bias",
        "eventID",
        "iMon",
        "iSys",
        "sciNum",
        "telNum",
        "tempSipm",
        "timestamp",
        "timestampEvt",
        "utc",
    ]
    for key in keys:
        if key in sci.keys():
            np.save(os.path.join(dir, f"{key}_{name}.npy"), sci[key])
        elif key in tel.keys():
            np.save(os.path.join(dir, f"{key}_{name}.npy"), tel[key])
        else:
            print(f"key {key} not found")


def data_load(path):
    sci = {}.fromkeys(["amp", "timestampEvt", "eventID", "sciNum"])
    tel = {}.fromkeys(
        ["bias", "iMon", "iSys", "telNum", "tempSipm", "timestamp", "utc"]
    )
    files = os.listdir(path)
    for file in files:
        header = file.split("_")[0]
        if header in sci.keys():
            sci[header] = np.load(os.path.join(path, file))
        elif header in tel.keys():
            tel[header] = np.load(os.path.join(path, file))
        elif header == "tempSiPM":
            tel["tempSipm"] = np.load(os.path.join(path, file))
        else:
            print(f"key {header} not found")
    return sci, tel
