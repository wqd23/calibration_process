import lib_reader as ver
from . import util_lib as util
from lib_reader.reader05.my_type import *
from dataclasses import dataclass, field
from .fitting import peak_fit
import numpy as np


def config_import(type, config):
    if isinstance(config, dict):
        return type(**config)
    elif isinstance(config, type):
        return config
    else:
        raise TypeError(
            f"config should be dict or {type.__name__}, but got {type(config).__name__}"
        )


@dataclass
class Read_config:
    # path of input file, '' for no bkg
    path: str = ""
    # ending param of dataReadout, for different file; 05B: 'normal' or 'xray'; 03B:'03b'
    ending: str = "normal"
    # configfile param of dataReadout, '' for search in the same dir
    config_file: str = ""
    # time cut
    time_cut: Optional[List[List[float]]] = None
    # other param for read_out
    kwarg: dict = field(default_factory=dict)


@dataclass
class Spectrum_config:
    # function to get temp bias corr factor
    corr: List[Callable[[float, float], float]] = field(
        default_factory=lambda: [lambda t, b: 1] * 4
    )
    # how to do temp bias corr, 'amp' for do corr for amp and get spectrum, 'spectrum' for get spectrum first and then do corr
    corr_style: str = "amp"
    adc_max: float = 16384.0
    bin_width: int = 2
    # how to change count into count rate,'' for simply divide (time[-1]-time[0])
    rate_style: str = ""
    kwarg: dict = field(default_factory=dict)


@dataclass
class Fit_config:
    fit_range: List[List[float]] = field(default_factory=lambda: [[None, None]] * 4)
    bkg_form: str = "lin"


class File_operation_05b:
    """settings and functions for transfer 4 channel data file into data point with error"""

    def __init__(
        self,
        path: str,
        read_config: Union[Read_config, Dict[str, Any]] = {},
        bkg_read_config: Union[Read_config, Dict[str, Any]] = {},
        spectrum_config: Union[Spectrum_config, Dict[str, Any]] = {},
        fit_config: Union[Fit_config, Dict[str, Any]] = {},
        nocache=False,
    ) -> None:
        self.path = path
        # config import
        self.read_config = config_import(Read_config, read_config)
        self.read_config.path = path
        self.bkg_read_config = config_import(Read_config, bkg_read_config)
        self.spectrum_config = config_import(Spectrum_config, spectrum_config)
        self.fit_config = config_import(Fit_config, fit_config)
        # read_config use a mutable {} as default value, read_config should never be changed
        del read_config, bkg_read_config, spectrum_config, fit_config

        # data readout=
        self.sci, self.tel = self.read_out(nocache=nocache)
        # # get spectrum
        # self.get_spectrum()
        # # peak fit
        # self.peak_fit()

    def __time_cut(
        self,
        data: Tuple[Dict[str, Float_array_4channel], Dict[str, Float_array_4channel]],
        time_cut: List[List[float]],
    ) -> Tuple[Dict[str, Float_array_4channel], Dict[str, Float_array_4channel]]:
        sci, tel = data
        for i in range(4):
            q1 = (sci["timestampEvt"][i] > time_cut[i][0]) * (
                sci["timestampEvt"][i] < time_cut[i][1]
            )
            for key in sci.keys():
                sci[key][i] = sci[key][i][q1]
            q2 = (tel["timestamp"][i] > time_cut[i][0]) * (
                tel["timestamp"][i] < time_cut[i][1]
            )
            if not np.any(q2):
                # if q2 is empty, maybe tel['timestamp'][i] is all 0.
                q2 = (tel["timestamp"][i] >= 0) * (tel["timestamp"][i] < time_cut[i][1])
            for key in tel.keys():
                tel[key][i] = tel[key][i][q2]
        return sci, tel

    def __read(
        self, config: Read_config, nocache: bool = False
    ) -> Tuple[Dict[str, Float_array_4channel], Dict[str, Float_array_4channel]]:
        if config.ending == "normal":
            data = ver.single_read05b_normal(config.path, overwrite_cache=nocache)
        elif config.ending == "xray":
            data = ver.single_read05b_xray(
                config.path, config.config_file, overwrite_cache=nocache
            )
        elif config.ending == "03b":
            data = ver.single_read03b(
                config.path, config.config_file, overwrite_cache=nocache
            )
        elif config.ending == "03b-src":
            data = ver.src_read03b(
                config.path, config.config_file, overwrite_cache=nocache
            )
        elif config.ending == "07":
            data = ver.single_read07(config.path, overwrite_cache=nocache)
        elif config.ending == "04":
            data = ver.single_read04(config.path, overwrite_cache=nocache)
        elif config.ending == "10b":
            data = ver.single_read10(config.path, overwrite_cache=nocache)
        elif config.ending == "11b":
            data = ver.single_read11(
                config.path, config.kwarg.get("mode", "wf"), overwrite_cache=nocache
            )
        else:
            raise ValueError(f"ending {config.ending} not supported")
        if config.time_cut != None:
            data = self.__time_cut(data, config.time_cut)
        return data

    def read_out(self, nocache=False):
        if self.bkg_read_config.path:
            self.bkg_sci, self.bkg_tel = self.__read(
                self.bkg_read_config, nocache=nocache
            )
        return self.__read(self.read_config, nocache=nocache)

    def __raw_spectrum(self, amp, bin_width, spec_range):
        spectrum, spectrum_err, x = util.count(
            amp, bin_width=bin_width, spec_range=spec_range
        )
        spectrum = [s / bin_width for s in spectrum]
        spectrum_err = [s / bin_width for s in spectrum_err]
        return spectrum, spectrum_err, x

    def __raw_rate(self, spectrum, spectrum_err, timestampEvt):
        if self.read_config.time_cut == None:
            time_spec = np.array(
                [np.max(np.hstack(timestampEvt)) - np.min(np.hstack(timestampEvt))] * 4
            )
        else:
            time_spec = np.array([np.max(time) - np.min(time) for time in timestampEvt])
        self.time = time_spec
        spectrum = [s / t for s, t in zip(spectrum, time_spec)]
        spectrum_err = [s / t for s, t in zip(spectrum_err, time_spec)]
        return spectrum, spectrum_err, 1 / time_spec, 1 / time_spec

    def get_spectrum(self):
        # corr = util.tb_corr(self.spectrum_config.corr, self.tel['tempSipm'], self.tel['bias'])
        # corr, corr_err = basic.tempBiasCorrection(self.tel['tempSipm'], self.tel['bias'], doCorr=True, ver = '03')
        corr = [
            c_f(np.mean(temp), np.mean(bias))
            for c_f, temp, bias in zip(
                self.spectrum_config.corr, self.tel["tempSipm"], self.tel["bias"]
            )
        ]
        amp = [factor * a for factor, a in zip(corr, self.sci["amp"])]
        # Integral of the spectrum is total count
        spectrum, spectrum_err, x = self.__raw_spectrum(
            amp,
            self.spectrum_config.bin_width,
            spec_range=[[0, c * self.spectrum_config.adc_max] for c in corr],
        )
        # integral is total count rate
        spectrum, spectrum_err, rate, rate_err = self.__raw_rate(
            spectrum, spectrum_err, self.sci["timestampEvt"]
        )

        if self.bkg_read_config.path:
            bkg_corr = util.tb_corr(
                self.spectrum_config.corr,
                self.bkg_tel["tempSipm"],
                self.bkg_tel["bias"],
            )
            bkg_amp = [factor * a for factor, a in zip(bkg_corr, self.bkg_sci["amp"])]
            # keep the same spec_range with data
            bkg_spectrum, bkg_spectrum_err, _ = self.__raw_spectrum(
                bkg_amp,
                self.spectrum_config.bin_width,
                spec_range=[[0, c * self.spectrum_config.adc_max] for c in corr],
            )
            bkg_spectrum, bkg_spectrum_err, bkg_rate, bkg_rate_err = self.__raw_rate(
                bkg_spectrum, bkg_spectrum_err, self.bkg_sci["timestampEvt"]
            )

            # get pure spectrum without bkg
            spectrum = [s - b for s, b in zip(spectrum, bkg_spectrum)]
            spectrum_err = [
                np.sqrt(s_err**2 + (b_s * b_rate_err) ** 2 + (b_err * b_rate) ** 2)
                for s_err, b_s, b_err, b_rate, b_rate_err in zip(
                    spectrum_err, bkg_spectrum, bkg_spectrum_err, bkg_rate, bkg_rate_err
                )
            ]

        self.spectrum = spectrum
        self.spectrum_err = spectrum_err
        self.x = x

    def peak_fit(self):
        fit_result = []
        for ich, (x, spectrum, spectrum_err, x_range, time) in enumerate(
            zip(
                self.x,
                self.spectrum,
                self.spectrum_err,
                self.fit_config.fit_range,
                self.time,
            )
        ):
            if x_range == None:
                fit_result.append(None)
                continue
            q = (x >= x_range[0]) * (x <= x_range[1])
            try:
                result = util.peak_fit(
                    x[q],
                    spectrum[q],
                    spectrum_err[q],
                    bkg_form=self.fit_config.bkg_form,
                )
            except util.FitError as e:
                print(e.args[0].fit_report())
                print(e.args[1])
                e.args[0].plot()
                raise util.FitError(f"failed to fit {self.path} channel {ich}")
            rate = result["peak_amplitude"]
            rate_err = np.sqrt(rate / time)
            fit_result.append(
                {
                    "a": result["peak_amplitude"],
                    "b": result["peak_center"],
                    "c": result["peak_sigma"],
                    "a_err": result["peak_amplitude_err"],
                    "b_err": result["peak_center_err"],
                    "c_err": result["peak_sigma_err"],
                    "rate": rate,
                    "rate_err": rate_err,
                    "resolution": 2
                    * np.sqrt(2 * np.log(2))
                    * result["peak_sigma"]
                    / result["peak_center"],
                    "resolution_err": 2
                    * np.sqrt(2 * np.log(2))
                    * np.sqrt(
                        (result["peak_sigma_err"] / result["peak_center"]) ** 2
                        + (
                            result["peak_sigma"]
                            * result["peak_center_err"]
                            / result["peak_center"] ** 2
                        )
                        ** 2
                    ),
                    "bkg": result["bkg"],
                }
            )
        self.fit_result = fit_result

    def save(self, path):
        # save after fit
        data = {
            "file": self.path,
            "fit_result": self.fit_result,
            "spectrum": self.spectrum,
            "x": self.x,
            "tel": self.tel,
            "config": {
                "read_config": self.read_config,
                "spectrum_config": self.spectrum_config,
                "fit_config": self.fit_config,
            },
        }
        util.pickle_save(data, path)


@dataclass
class Fit_cfg4ch:
    """settings for 4 channel fit"""

    fit_range: List[List[float]] = field(default_factory=lambda: [[None, None]] * 4)
    bkg_form: List[str] = field(default_factory=lambda: ["lin"] * 4)
    peak_form: List[str] = field(default_factory=lambda: ["gaus"] * 4)


class File_operation_10b(File_operation_05b):
    """settings and functions for transfer 4 channel data file into data point with error, with compton fitting ability"""

    def __init__(
        self,
        read_config: Union[Read_config, Dict[str, Any]] = {},
        bkg_read_config: Union[Read_config, Dict[str, Any]] = {},
        spectrum_config: Union[Spectrum_config, Dict[str, Any]] = {},
        fit_config: Union[Fit_cfg4ch, Dict[str, Any]] = {},
    ):
        path = (
            read_config.path
            if isinstance(read_config, Read_config)
            else read_config["path"]
        )
        super().__init__(path, read_config, bkg_read_config, spectrum_config)
        self.fit_config = config_import(Fit_cfg4ch, fit_config)
        # read_config use a mutable {} as default value, read_config should never be changed
        del read_config, bkg_read_config, spectrum_config, fit_config
        self.sci, self.tel = self.read_out()

    def peak_fit(self, **kwargs):
        fit_result = []
        for ich, (x, spectrum, spectrum_err, x_range, time, bkg, peak) in enumerate(
            zip(
                self.x,
                self.spectrum,
                self.spectrum_err,
                self.fit_config.fit_range,
                self.time,
                self.fit_config.bkg_form,
                self.fit_config.peak_form,
            )
        ):
            if x_range == None:
                fit_result.append(None)
                continue
            q = (x >= x_range[0]) * (x <= x_range[1])

            try:
                result = peak_fit(
                    x[q],
                    spectrum[q],
                    spectrum_err[q],
                    peak_form=peak,
                    bkg_form=bkg,
                    **kwargs,
                )
            except Exception as e:
                print(e)
                raise util.FitError(f"failed to fit {self.path} channel {ich}")
            rate = result["peak_amplitude"]
            rate_err = np.sqrt(rate / time)
            fit_result.append(
                {
                    "a": result["peak_amplitude"],
                    "b": result["peak_center"],
                    "c": result["peak_sigma"],
                    "a_err": result["peak_amplitude_err"],
                    "b_err": result["peak_center_err"],
                    "c_err": result["peak_sigma_err"],
                    "rate": rate,
                    "rate_err": rate_err,
                    "resolution": result["peak_resolution"],
                    "resolution_err": result["peak_resolution_err"],
                    "bkg": result["bkg"],
                }
            )
        self.fit_result = fit_result
