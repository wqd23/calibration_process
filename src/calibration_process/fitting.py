from lmfit import Model
from lmfit.models import fwhm_expr, height_expr, COMMON_INIT_DOC, COMMON_GUESS_DOC, GaussianModel, LinearModel, ExponentialModel, QuadraticModel
from scipy.special import erfc
import numpy as np

def get_peak_mod(name:str, *args, **kwargs):
    peak_mod = {
        "gaus": GaussianModel,
        "comp": ComptonModelSimple
    }
    return peak_mod.get(name)(*args, **kwargs)


def get_bkg_mod(name:str, *args, **kwargs):
    bkg_mod = {
        "gaus": GaussianModel,
        "lin": LinearModel,
        "exp": ExponentialModel,
        "quad": QuadraticModel
    }
    return bkg_mod.get(name)(*args, **kwargs)




def compton_erfc(x, center=1700, sigma=100, amplitude=0.5e-9, b=0):
    e = (x-center)/sigma

    return amplitude*erfc(e) + b


class ComptonModelSimple(Model):
    r"""A simple model based on erfc function to fit compton edge.

    """

    fwhm_factor = 2*np.sqrt(2*np.log(2))
    height_factor = 1./np.sqrt(2*np.pi)

    def __init__(self, independent_vars=['x'], prefix='', nan_policy='raise',
                 **kwargs):
        kwargs.update({'prefix': prefix, 'nan_policy': nan_policy,
                       'independent_vars': independent_vars})
        super().__init__(compton_erfc, **kwargs)
        self._set_paramhints_prefix()

    def _set_paramhints_prefix(self):
        self.set_param_hint('sigma', min=0)
        self.set_param_hint('fwhm', expr=fwhm_expr(self))
        self.set_param_hint('height', expr=height_expr(self))


    def guess(self, data, x, **kwargs):
        """Estimate initial model parameter values from data."""
        center = (x[0] + x[-1]) / 2
        amplitude = (data[0] + data[-1]) / 2
        sigma = (x[-1] - x[0]) / 4
        pars = self.make_params(center=center, amplitude=amplitude, sigma=sigma, b = 0)
        pars[f'{self.prefix}sigma'].set(min=0.0)
        return pars
    __init__.__doc__ = COMMON_INIT_DOC
    guess.__doc__ = COMMON_GUESS_DOC