from lmfit import Model
from lmfit.models import (
    fwhm_expr,
    height_expr,
    COMMON_INIT_DOC,
    COMMON_GUESS_DOC,
    GaussianModel,
    LinearModel,
    ExponentialModel,
    QuadraticModel,
)
from scipy.special import erfc
from lib_reader.reader05.my_type import *

import numpy as np


def get_peak_mod(name: str, *args, **kwargs) -> Model:
    peak_mod = {"gaus": GaussianModel, "compErfc": ComptonModelSimple}
    mod = peak_mod.get(name)
    if mod is None:
        raise ValueError(
            f"Unknown peak model name: {name}. Available options are: {list(peak_mod.keys())}"
        )
    return mod(*args, **kwargs)


def get_bkg_mod(name: str, *args, **kwargs) -> Model:
    bkg_mod = {
        "gaus": GaussianModel,
        "lin": LinearModel,
        "exp": ExponentialModel,
        "quad": QuadraticModel,
    }
    mod = bkg_mod.get(name)
    if mod is None:
        raise ValueError(
            f"Unknown background model name: {name}. Available options are: {list(bkg_mod.keys())}"
        )

    return mod(*args, **kwargs)


def compton_erfc(x, center=1700, sigma=100, amplitude=0.5e-9, b=0):
    e = (x - center) / sigma

    return amplitude * erfc(e) + b


class ComptonModelSimple(Model):
    r"""A simple model based on erfc function to fit compton edge."""

    fwhm_factor = 2 * np.sqrt(2 * np.log(2))
    height_factor = 1.0 / np.sqrt(2 * np.pi)

    def __init__(self, independent_vars=["x"], prefix="", nan_policy="raise", **kwargs):
        kwargs.update(
            {
                "prefix": prefix,
                "nan_policy": nan_policy,
                "independent_vars": independent_vars,
            }
        )
        super().__init__(compton_erfc, **kwargs)
        self._set_paramhints_prefix()

    def _set_paramhints_prefix(self):
        self.set_param_hint("sigma", min=0)
        self.set_param_hint("fwhm", expr=fwhm_expr(self))
        self.set_param_hint("height", expr=height_expr(self))

    def guess(self, data, x, **kwargs):
        """Estimate initial model parameter values from data."""
        center = (x[0] + x[-1]) / 2
        amplitude = 2 * data[0]
        sigma = (x[-1] - x[0]) / 4
        pars = self.make_params(center=center, amplitude=amplitude, sigma=sigma, b=0)
        pars[f"{self.prefix}sigma"].set(min=0.0)
        return pars

    __init__.__doc__ = COMMON_INIT_DOC
    guess.__doc__ = COMMON_GUESS_DOC


def peak_fit(
    x: Float1D,
    data: Float1D,
    error: Float1D,
    peak_form: str,
    bkg_form: Optional[str] = None,
    **kwargs,
) -> Dict[str, Any]:
    peak_mod = get_peak_mod(peak_form, prefix="peak_")
    param = peak_mod.guess(data, x)
    mod = peak_mod
    if bkg_form is not None:
        bkg_mod = get_bkg_mod(bkg_form, prefix="bk_")
        mod += bkg_mod
        param.update(bkg_mod.make_params())
    param["peak_amplitude"].min = 0
    param["peak_center"].min, param["peak_center"].max = x[0], x[-1]
    param["peak_sigma"].min, param["peak_sigma"].max = 0, (x[-1] - x[0]) / 3.0
    param.add("peak_resolution", expr="peak_fwhm/peak_center")
    for k, v in kwargs.items():
        if k in param:
            param[k].set(value=v)
    result = mod.fit(data, param, x=x, weights=1.0 / error)
    fit_result = {k: v for k, v in result.values.items() if k.startswith("peak_")}
    fit_result.update(
        {
            "peak_amplitude_err": result.params["peak_amplitude"].stderr,
            "peak_center_err": result.params["peak_center"].stderr,
            "peak_sigma_err": result.params["peak_sigma"].stderr,
            "peak_resolution_err": result.params["peak_resolution"].stderr,
        }
    )
    if bkg_form is None:
        bkg = {"bkg_info": "None", "bkg_func": lambda x: 0 * x}
    else:
        bkg = {
            "bkg_info": bkg_form,
            "bkg_func": lambda x: bkg_mod.eval(x=x, **result.best_values),
        }
        bkg.update({k: v for k, v in result.values.items() if k.startswith("bk_")})
    fit_result["bkg"] = bkg
    fit_result["peak_func"] = lambda x: peak_mod.eval(x=x, **result.best_values)
    return fit_result
