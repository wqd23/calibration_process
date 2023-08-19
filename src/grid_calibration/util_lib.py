# -*- coding:utf-8 -*-
"""
GRID version independent Functions Library
----------

Basic functions for GRID data processing\n
"""
import lib_reader.reader05.gridBasicFunctions as basic
from lib_reader.reader05.my_type import *
import nptyping as npt

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import lmfit.models as fitmodel
import lmfit
from scipy.optimize import curve_fit
import dill as pickle
import json
import os
import datetime
import math
def count_spectrum(amp, nbins, spec_range, bin_width, adc_max):
    spectrum, x = basic.getSpectrum(
        amp, nbins=nbins, specRange=spec_range, binWidth=bin_width, adcMax=adc_max)
    spectrum_err: Float_array_4channel = [basic.gehrelsErr(
        spectrum_ich) for spectrum_ich in spectrum]  # type: ignore
    return spectrum, spectrum_err, x
def count(amp:Float_array_4channel, bin_width:float, spec_range:List[List[float]])->Tuple[Float_array_4channel, Float_array_4channel, Float_array_4channel]:
    spectrum_all = []
    spectrum_err_all = []
    x_all = []
    for a,s_range in zip(amp, spec_range):
        binEdges = np.arange(s_range[0], s_range[1] + bin_width, bin_width)
        spectrum, x = np.histogram(a, bins = binEdges)
        assert spectrum.sum() == a.size, "spectrum sum is not equal to the number of events"
        x = (x[:-1] + x[1:]) / 2
        spectrum_err = basic.gehrelsErr(spectrum)
        spectrum_all.append(spectrum)
        spectrum_err_all.append(spectrum_err)
        x_all.append(x)
    return spectrum_all, spectrum_err_all, x_all

def raw_plot_single_channel(spectrum, x, title, ax=None, **kwargs):
    if ax is None:
        ax = plt.subplot()
    ax.step(x, spectrum, where='mid', label="raw data")
    ax.set_xlabel('ADC/channel')
    ax.set_ylabel('count rate/cps')
    ax.grid()
    ax.set_title(title)
    # kwargs
    # log scale
    if kwargs.get('log_scale', 'linear') in ["log", "semilogx"]:
        ax.set_xscale('log')
    if kwargs.get('log_scale', 'linear') in ["log", "semilogy"]:
        ax.set_yscale('log')
    return ax

def __raw_plot(spectrum:Float_array_4channel, x:Float_array_4channel,title:str, **kwargs):
    default_kwargs = {"log_scale": "linear"}
    kwargs = {**default_kwargs, **kwargs}
    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(4, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
    for ich, (spectrum_ich, x_ich) in enumerate(zip(spectrum, x)):
        ax = fig.add_subplot(gs[ich])
        if ich == 0:
            ax.set_title(title)
        # spectrum plot
        plt.step(x_ich, spectrum_ich,where='mid',label="raw data")
        ax.set_xlabel('ADC/channel')
        ax.set_ylabel('count rate/cps')
        ax.grid()
        # ax.set_xlim((0, 1000))
        # kwargs
        # log scale
        if kwargs["log_scale"] in ["log", "semilogx"]:
            ax.set_xscale('log')
        if kwargs["log_scale"] in ["log", "semilogy"]:
            ax.set_yscale('log')
        if kwargs.get("x_lim", None) is not None:
            ax.set_xlim(kwargs['x_lim'])
    return fig, fig.axes
def raw_plot(spectrum:Float_array_4channel, x:Float_array_4channel,title:str,save_path: Optional[str] = None, **kwargs):
    fig, axs = __raw_plot(spectrum, x, title, **kwargs)
    if save_path != None:
        fig.savefig(save_path)
        plt.close(fig)
    else:
        plt.show()

def fit_plot_single_channel(spectrum, x, fit_result, title, bkgForm, fit_range=None,ax=None, **kwargs):
    if ax is None:
        ax = plt.subplot()
    if fit_result is None:
        ax.text(x[np.where(spectrum==np.max(spectrum))],max(spectrum),"None fit range offered\n no fit for this channel")
        return ax
    raw_plot_single_channel(spectrum, x, title, ax=ax, **kwargs)
    amplitude = fit_result['a']
    center = fit_result['b']
    sigma = fit_result['c']
    fit_peak = basic.gaussianFunction([amplitude, center, sigma], x)
    bkg = fit_result['bkg']
    fit_bkg = bkg['bkg_func'](x)
    fit_total = fit_peak+fit_bkg
    if fit_range is None:
        plot_range = [center-3*sigma,center+3*sigma]
    else:
        plot_range = fit_range
    qplot = (x >= plot_range[0]) * (x < plot_range[1])
    ax.plot(x[qplot], fit_peak[qplot],label = 'Gaussian Peak Fit')
    ax.plot(x[qplot], fit_bkg[qplot],label = 'Background Fit')
    ax.plot(x[qplot], fit_total[qplot],label = 'Gaussian Peak and Background Fit')
    info = f"""
    count rate = {fit_result["rate"]:.2e} $\\pm${fit_result["rate_err"]:.3e}
    center = {fit_result["b"]:.3e} $\\pm${fit_result["b_err"]:.3e}
    resolution = {fit_result["resolution"]:.3e} $\\pm${fit_result["resolution_err"]:.3e}
    """
    plt.text(center - 0.5 * sigma, max(fit_peak[qplot]) * 0.25,  info,
            fontsize = 10, bbox = dict(facecolor = 'pink', alpha = 0.1), horizontalalignment = 'center', verticalalignment = 'center')
    ax.legend(loc = 0)
    if fit_range is not None:
        ax.axvline(plot_range[0])
        ax.axvline(plot_range[1])
    # set x,y lim
    ax.set_xlim(plot_range[0]-3*sigma, plot_range[1]+3*sigma)
    ax.set_ylim(0, 1.2 * np.max(spectrum[qplot]))
    return ax
def __fit_plot(spectrum:Float_array_4channel, x:Float_array_4channel, fit_result:List[Union[Dict[str, float],None]], title:str, bkgForm:str, fit_range = None, **kwargs):        
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
    
    bkg_func = lambda param, x: (param[0] * x + param[1]) if bkgForm == 'lin' else basic.quadFunction(param, x)

    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(4, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
    for ich, (spectrum_ich, x_ich, fit_ich) in enumerate(zip(spectrum, x, fit_result)):
        ax = fig.add_subplot(gs[ich])
        if ich == 0:
            ax.set_title(title)
        # spectrum plot
        plt.step(x_ich, spectrum_ich,where='mid',label="raw data")
        if fit_ich != None:
            # fit plot data
            amplitude = fit_ich['a']
            center = fit_ich['b']
            sigma = fit_ich['c']
            fit_peak = basic.gaussianFunction([amplitude, center, sigma], x_ich)
            bkg = fit_ich['bkg']
            fit_bkg = bkg['bkg_func'](x_ich)
            fit_total = fit_peak+fit_bkg
            # fit plot
            if fit_range != None and fit_range[ich] != None:
                # fit range set
                plot_range =  [fit_range[ich][0], fit_range[ich][1]]
            else:
                # fit range not set
                plot_range =  [center-3*sigma,center+3*sigma]
            qplot = (x_ich >= plot_range[0]) * (x_ich < plot_range[1])
            plt.plot(x_ich[qplot], fit_peak[qplot],label = 'Gaussian Peak Fit')
            plt.plot(x_ich[qplot], fit_bkg[qplot],label = 'Background Fit')
            plt.plot(x_ich[qplot], fit_total[qplot],label = 'Gaussian Peak and Background Fit')
            # fit info
            info = f"""
            count rate = {fit_ich["rate"]:.2e} $\\pm${fit_ich["rate_err"]:.3e}
            center = {fit_ich["b"]:.3e} $\\pm${fit_ich["b_err"]:.3e}
            resolution = {fit_ich["resolution"]:.3e} $\\pm${fit_ich["resolution_err"]:.3e}
            """
            plt.text(center - 0.5 * sigma, max(fit_peak[qplot]) * 0.25,  info,
                    fontsize = 10, bbox = dict(facecolor = 'pink', alpha = 0.1), horizontalalignment = 'center', verticalalignment = 'center')
            # set x,y lim
            ax.set_xlim(plot_range[0]-3*sigma, plot_range[1]+3*sigma)
            ax.set_ylim(0, 1.2 * np.max(spectrum_ich[qplot]))
        else:
            plt.text(x_ich[np.where(spectrum_ich==np.max(spectrum_ich))],max(spectrum_ich),"None fit range offered\n no fit for this channel")

        ax.set_xlabel('ADC/channel')
        ax.set_ylabel('count rante/cps')
        ax.legend(loc = 0)
        ax.grid()
        # fit range
        if fit_range != None and fit_range[ich] != None:
            plt.axvline(fit_range[ich][0])
            plt.axvline(fit_range[ich][1])
        # kwargs
        # log scale
        if kwargs["log_scale"] in ["log", "semilogx"]:
            ax.set_xscale('log')
        if kwargs["log_scale"] in ["log", "semilogy"]:
            ax.set_yscale('log')
    return fig, fig.axes
def fit_plot(spectrum:Float_array_4channel, x:Float_array_4channel, fit_result:List[Union[Dict[str, float],None]], title:str, bkgForm:str, fit_range = None,save_path: Optional[str] = None, **kwargs):
    fig, axs = __fit_plot(spectrum, x, fit_result, title, bkgForm, fit_range, **kwargs)
    if save_path != None:
        fig.savefig(save_path)
        plt.close(fig)
    else:
        plt.show()

class FitError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
        # message for print
        self.message = args[-1]


def peak_fit(x:Float1D, data:Float1D, error:Float1D, bkg_form:Optional[str]=None)-> Dict[str, Any]:
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
    mod = fitmodel.GaussianModel(prefix='peak_')
    param = mod.guess(data, x = x)
    if bkg_form == 'lin':
        b_mod = lmfit.models.LinearModel(prefix = 'bk_')
    elif bkg_form == 'quad':
        b_mod = lmfit.models.QuadraticModel(prefix = 'bk_')
    elif bkg_form == 'exp':
        b_mod = lmfit.models.ExponentialModel(prefix = 'bk_')
    else:
        b_mod = fitmodel.GaussianModel(prefix='bk_')
    if bkg_form != None:
        # while use bkg, right half to guess data
        data_param = mod.guess(data[len(data)//2:], x = x[len(data)//2:])
        bkg_data = data - mod.eval(data_param, x=x)
        b_param = b_mod.guess(bkg_data, x=x)
        mod = mod + b_mod
        param = data_param + b_param
    # set bound for param
    # assume center is in the xdata range, xdata range bigger than 3*sigma(half the peak)
    # param['peak_amplitude'].min, param['peak_amplitude'].max = 0, 10*np.max(ydata)
    param['peak_amplitude'].min = 0
    param['peak_center'].min, param['peak_center'].max = x[0], x[-1]
    param['peak_sigma'].min, param['peak_sigma'].max = 0, (x[-1]-x[0])/3.
    param.add('peak_resolution', expr='peak_fwhm/peak_center')
    result = mod.fit(data,param, x=x,weights=1./error)
    fitResult = {k:v for k,v in result.values.items() if 'bk_' not in k}
    fitResult.update({
            'peak_amplitude_err':           result.params['peak_amplitude'].stderr,
            'peak_center_err':                result.params['peak_center'].stderr,
            'peak_sigma_err':                 result.params['peak_sigma'].stderr,
    })

    if bkg_form == 'lin':
        fitResult.update({'bk_a':       result.params['bk_slope'].value, 
                            'bk_b':        result.params['bk_intercept'].value, 
            })
    if bkg_form != None:
        fitResult['bkg'] = {
            'bkg_info': bkg_form,
            'bkg_func': lambda x: b_mod.eval(x=x, **result.best_values),
        }
        fitResult['bkg'].update({k:v for k,v in result.values.items() if 'bk_' in k})
    else:
        fitResult['bkg'] = {
            'bkg_info': 'None',
            'bkg_func': lambda x: 0*x,   
        }
    # fit failed or can't get uncertainty
    if not result.success or any(map(lambda x: x== None, fitResult.values())):
        # param too much, that can't get uncertainty
        if bkg_form == 'lin':
            # fit without bkg
            try:
                fitResult = peak_fit(x, data, error)
            except FitError as e:
                raise FitError(e.args[0], f"failed to fit with {bkg_form} background and failed to fit without background")
            # add the linear zero bkg
            fitResult.update({'bk_a':       0, 
                            'bk_b':        0, 
            })
        else:
            raise FitError(result,f"failed to fit with {bkg_form} background")
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
    xdata = x[:,0]
    ydata = x[:,1]
    Vov = ydata - k * xdata - V0
    return G0 * Vov ** 2 * (-xdata ** 2 + b * xdata + c)
def tempbias2DFunction(temp, bias, G0, k, V0, b, c):
    Vov = bias - k * temp - V0
    return G0 * Vov ** 2 * (-temp ** 2 + b * temp + c)
def resolutionFunction(x, a,b,c):
    y = np.sqrt(a * x * x + b * x + c) / x
    y[a * x * x + b * x + c < 0] = 0.
    return y
def temp_bias_lmfit(center:Float1D, center_err:Float1D, temp:Float1D, temp_err:Float1D,bias:Float1D,bias_err:Float1D):

    data = np.stack([temp,bias, center], axis=1)
    error = np.stack([temp_err,bias_err, center_err], axis=1)
    model = lmfit.Model(tempbias2DFunctionInternal)
    # param = model.guess(data[:,2], data[:,0:2])

    model.set_param_hint('G0', value=0.030317363114139174,      min=1e-5,  max=100.)
    model.set_param_hint('k', value=18.9e-3, min=1e-3,  max=50e-3)
    model.set_param_hint('V0', value=24.6,   min=23.,   max=25.)
    model.set_param_hint('b', value=0.,      min=-500., max=500.)
    model.set_param_hint('c', value=1e4,     min=1e3,   max=1e8)
    param = model.make_params()
    result = model.fit(data[:,2], x = data[:,0:2], params=param, weights=1./error[:,2], method='trf')
    fitResult = {
                'G0':             result.best_values['G0'], 
                'k':               result.best_values['k'], 
                'V0':             result.best_values['V0'], 
                'b':               result.best_values['b'], 
                'c':               result.best_values['c'], 
                'G0_err':             param['G0'].stderr, 
                'k_err':               param['k'].stderr, 
                'V0_err':             param['V0'].stderr, 
                'b_err':               param['b'].stderr, 
                'c_err':               param['c'].stderr, 
                'chisquare':          result.chisqr
    }
    if not result.success or any(map(lambda x: x== None, fitResult.values())):
        raise FitError(result, "faled to do temp bias 2d fit")
def temp_bias_fit_curvefit(center:Float1D, center_err:Float1D, temp:Float1D,bias:Float1D)-> Dict[str, float]:
    data = np.stack([temp, bias], axis=1)
    initial_guess = [0.00078,18.9e-3,24.6,0,1e4]
    try:
        popt, pcov = curve_fit(tempbias2DFunctionInternal, data, center, sigma = center_err, absolute_sigma = True, p0=initial_guess)
    except RuntimeError as e:
        raise FitError(e.args[0])
    perr = np.sqrt(np.diag(pcov))
    param = popt.tolist()
    chisq = sum((basic.residualTempbias2D(param, temp, bias, center) / center_err) ** 2) / (len(center_err) - len(param))
    fitResult = {
                    'G0':             param[0], 
                    'k':               param[1], 
                    'V0':             param[2], 
                    'b':               param[3], 
                    'c':               param[4], 
                    'G0_err':             perr[0], 
                    'k_err':               perr[1], 
                    'V0_err':             perr[2], 
                    'b_err':               perr[3], 
                    'c_err':               perr[4], 
                    'chisquare':          chisq, 
        }
    return fitResult

def plot_3d(xdata, ydata, zdata, xx,yy,zz, label, ax, title=None):
    ax.scatter3D(xdata, ydata, zdata, label="data point") # type: ignore
    ax.plot_surface(xx,yy,zz, label="fit", alpha=0.4) # type: ignore

    ax.set_title(title)
    ax.set_xlabel(label[0])
    ax.set_ylabel(label[1])
    ax.set_zlabel(label[2]) # type: ignore

def fit_err_plot_2d(data:npt.NDArray[npt.Shape["*,2"], npt.Floating], zdata:Float1D, func:Callable,label:Tuple[str,str,str], title:Optional[str] = None, save_path:Optional[str] = None):
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
    xdata = data[:,0]
    ydata = data[:,1]
    x = np.arange(np.min(xdata)-1, np.max(xdata)+1)
    y = np.arange(np.min(ydata)-1, np.max(ydata)+1)
    xx,yy = np.meshgrid(x,y)
    zz = func(np.stack((xx,yy), axis=1))
    fig_num = 2
    fig = plt.figure(figsize=(10, 5))
    if title != None:
        fig.suptitle(title)
    # 3d plot
    # ax = fig.add_subplot(fig_num,1,1,projection='3d')
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    plot_3d(xdata, ydata, zdata, xx,yy,zz,label, ax, title="3d plot")
    # 2d residual plot
    residual = func(data)-zdata
    ax = fig.add_subplot(1,2,2)
    cs = ax.tricontourf(xdata, ydata, residual/zdata*100) # type: ignore
    plt.colorbar(cs)
    ax.scatter(xdata, ydata, np.abs(residual/zdata*100), c='r', marker='o') # type: ignore
    ax.set_title("relative residual$\\times 100$")
    ax.set_xlabel(label[0])
    ax.set_ylabel(label[1])
    # set x,y lim to min-0.01(max-min), max+0.01(max-min), to make all data points in view
    ax.set_xlim(1.01*np.min(xdata)-0.01*np.max(xdata), 1.01*np.max(xdata)-0.01*np.min(xdata))
    ax.set_ylim(1.01*np.min(ydata)-0.01*np.max(ydata), 1.01*np.max(ydata)-0.01*np.min(ydata))
    plt.tight_layout()
    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)

def pickle_save(data, path:str):
    """save basic data with pickle

    Parameters
    ----------
    data : Any
        basic python data
    path : str
        full path to save the data
    """    
    with open(path, 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
def pickle_load(path)-> Any:
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
    with open(path, 'rb') as f:
        return pickle.load(f)

def json_save(data, path:str):
    """save the data into a json file

    Parameters
    ----------
    data : Any
        any oject could save as json
    path : str
        full path "path/to/file.json"
    """    
    with open(path, 'w') as f:
        f.write(json.dumps(data))

def json_time_save(data, path:str):
    """save data with a time appended file name, wrapper for json_save

    Parameters
    ----------
    data : Any
        data object
    path : str
        full path, time info will be appended to it
    """    
    json_save(data, f"{os.path.splitext(path)[0]}_{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}.json")

def json_headtime_save(data, path:str):
    """save data with a time appended file name, wrapper for json_save

    Parameters
    ----------
    data : Any
        data object
    path : str
        full path, time info will be appended to it
    """    
    json_save(data, f"{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}_{os.path.splitext(path)[0]}.json")

def time_save(path):

    return f"{os.path.splitext(path)[0]}_{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}{os.path.splitext(path)[1]}"
def headtime(path):
    return f"{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}_{path}"
def json_load(path:str):
    with open(path, 'r') as f:
        return json.loads(f.read())

def tb_corr(corr:List[Callable[[float, float],float]], temp:Float_array_4channel, bias:Float_array_4channel)->Float_4channel:
    temp_avg = [np.average(t) for t in temp]
    bias_avg = [np.average(b) for b in bias]
    corr_factor = [f(t,b) for f,t,b in zip(corr,temp_avg, bias_avg)]
    return corr_factor

def ec_plot(energy,center, result, src_energy, x_energy, src_result,x_result, save_path):
    plt.ion()
    q_low = energy<50.2
    q_high = energy>=50.2
    center_low = [c[q_low] for c in center]
    center_high = [c[q_high] for c in center]

    adc_low = [np.arange(np.min(center_low[i]),np.max(center_low[i])) for i in range(4)]
    adc_high = [np.arange(np.min(center_high[i]),np.max(center_high[i])) for i in range(4)]
    energy_low = [np.polyval(result[i]['EC_low'], adc_low[i]) for i in range(4)]
    energy_high = [np.polyval(result[i]['EC_high'], adc_high[i]) for i in range(4)]

    src_center = [np.array([fit[i]['b'] for fit in src_result]) for i in range(4)]
    src_center_err = [np.array([fit[i]['b_err'] for fit in src_result]) for i in range(4)]
    x_center = [np.array([fit[i]['b'] for fit in x_result]) for i in range(4)]
    x_center_err = [np.array([fit[i]['b_err'] for fit in x_result]) for i in range(4)]

    
    gs = gridspec.GridSpec(2, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95,height_ratios = [4, 1])
    for i in range(4):
        fig = plt.figure(figsize = (12, 8))
        ax = fig.add_subplot(gs[0])
        ax.errorbar(x_center[i], x_energy, xerr = x_center_err[i], fmt = 's', mfc = 'white', ms = 6, elinewidth = 1, capsize = 3, barsabove = True, zorder = 1, label = f' CH{i}')
        ax.errorbar(src_center[i], src_energy, xerr = src_center_err[i], fmt = '^', mfc = 'white', ms = 6, elinewidth = 1, capsize = 3, barsabove = True, zorder = 0, \
            label = f'source CH{i}')
        ax.plot(adc_low[i], energy_low[i], linestyle = '-', label = f'quadratic fit on EC data of ch{i}, < 50.2keV')
        ax.plot(adc_high[i], energy_high[i], linestyle = '-', label = f'quadratic fit on EC data of ch{i}, > 50.2keV')
        ax.set_xlabel('ADC')
        ax.set_ylabel('Energy/keV')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_xlim(80., 16384.0)
        ax.set_ylim(10., 1500.)
        ax.legend(loc = 0)
        ax.grid()
        fig.savefig(time_save(os.path.join(save_path, f"ec_fit_ch{i}.png")))
    
    e_low = np.arange(np.min(energy), 50.2)
    e_high = np.arange(50.2, np.max(energy))
    resolution_low = [resolutionFunction(e_low, result[i]['resolution_low'][0],result[i]['resolution_low'][1],result[i]['resolution_low'][2]) for i in range(4)]
    resolution_high = [resolutionFunction(e_high, result[i]['resolution_high'][0],result[i]['resolution_high'][1],result[i]['resolution_high'][2]) for i in range(4)]
    src_resolution = [np.array([fit[i]['resolution'] for fit in src_result]) for i in range(4)]
    src_resolution_err = [np.array([fit[i]['resolution_err'] for fit in src_result]) for i in range(4)]
    x_resolution = [np.array([fit[i]['resolution'] for fit in x_result]) for i in range(4)]
    x_resolution_err = [np.array([fit[i]['resolution_err'] for fit in x_result]) for i in range(4)]

    gs = gridspec.GridSpec(2, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95,height_ratios = [4, 1])
    for i in range(4):
        fig = plt.figure(figsize = (12, 8))
        ax = fig.add_subplot(gs[0])
        if not any(map(math.isinf,x_resolution_err[i])):
            ax.errorbar(x_energy, x_resolution[i]*100, yerr = x_resolution_err[i]*100, fmt = 's', mfc = 'white', ms = 6, elinewidth = 1, capsize = 3, barsabove = True, zorder = 1, label = f' CH{i}')
        else:
            ax.scatter(x_energy, x_resolution[i]*100,label = f' CH{i}, inf in error')
        if not any(map(math.isinf,src_resolution_err[i])):
            ax.errorbar(src_energy, src_resolution[i]*100, yerr = src_resolution_err[i]*100, fmt = '^', mfc = 'white', ms = 6, elinewidth = 1, capsize = 3, barsabove = True, zorder = 0, \
            label = f'source CH{i}')
        else:
            ax.scatter(src_energy, src_resolution[i]*100,label = f'source CH{i}, inf in error')
        ax.plot(e_low, resolution_low[i]*100, linestyle = '-', label = f'Fit on resolution data of ch{i}, < 50.2keV')
        ax.plot(e_high, resolution_high[i]*100, linestyle = '-', label = f'Fit on resolution data of ch{i}, > 50.2keV')
        ax.set_xlabel('Energy/keV')
        ax.set_ylabel('Resolution/%')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_xlim(20, 1500)
        ax.set_ylim(8, 150)
        ax.legend(loc = 0)
        ax.grid()
        fig.savefig(time_save(os.path.join(save_path, f"resolution_fit_ch{i}.png")))

def resolution_fit(energy:Float1D, resolution:Float1D, resolution_err:Float1D):
    # initial value
    data = (resolution*energy)**2
    error = 2*energy**2*resolution_err
    p0, pcov = np.polyfit(energy, data, deg=2, w=1./error, cov=True)
    perr = np.sqrt(np.diag(pcov))
    popt = list(p0)
    perr = list(perr)
    # mod = lmfit.Model(resolutionFunction)
    # param = mod.make_params(a=p0[0], b=p0[1], c=p0[2])
    # result = mod.fit(resolution, param, weights=1./resolution_err, x= energy)
    # popt = [result.best_values['a'], result.best_values['b'], result.best_values['c']]
    # perr = [param['a'].stderr, param['b'].stderr, param['c'].stderr]
    return popt, perr

def get_key(file_name:str):
    assert os.path.exists(file_name), f'{file_name} does not exist!'
    return os.path.basename(file_name)

def get_hpge_data(path):
    data = np.loadtxt(path)
    counts = data
    # unit: keV
    energy = np.arange(0, len(counts), 1)
    return energy, counts

def print_temp(tel):
    temp = tel['tempSipm']
    avg = [np.mean(temp[i]) for i in range(4)]
    std = [np.std(temp[i]) for i in range(4)]
    print(f'CH0: {avg[0]:.2f}+/-{std[0]:.2f} CH1: {avg[1]:.2f}+/-{std[1]:.2f} CH2: {avg[2]:.2f}+/-{std[2]:.2f} CH3: {avg[3]:.2f}+/-{std[3]:.2f}')

def data_save(sci, tel, path, name):
    dir = os.path.join(path, name)
    if not os.path.exists(dir):
        os.makedirs(dir)
    keys = ['amp', 'bias', 'eventID', 'iMon', 'iSys', 'sciNum', 'telNum', 'tempSipm', 'timestamp', 'timestampEvt', 'utc']
    for key in keys:
        if key in sci.keys():
            np.save(os.path.join(dir, f'{key}_{name}.npy'), sci[key])
        elif key in tel.keys():
            np.save(os.path.join(dir, f'{key}_{name}.npy'), tel[key])
        else:
            print(f'key {key} not found')

def data_load(path):
    sci = {}.fromkeys(['amp','timestampEvt', 'eventID','sciNum'])
    tel = {}.fromkeys(['bias', 'iMon', 'iSys', 'telNum', 'tempSipm', 'timestamp', 'utc'])
    files = os.listdir(path)
    for file in files:
        header = file.split('_')[0]
        if header in sci.keys():
            sci[header] = np.load(os.path.join(path, file))
        elif header in tel.keys():
            tel[header] = np.load(os.path.join(path, file))
        elif header == 'tempSiPM':
            tel['tempSipm'] = np.load(os.path.join(path, file))
        else:
            print(f'key {header} not found')
    return sci, tel
    