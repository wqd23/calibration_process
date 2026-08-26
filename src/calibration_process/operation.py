import lib_reader as ver
from . import util_lib as util
from lib_reader.reader05.my_type import *
from lib_plot import plot
import os
import numpy as np
from . import file_lib
from pathlib import Path


class Operation:
    """wrapper for tb_op, x of ec_op, and src of ec_op"""

    def __init__(self, files, file_config, op) -> None:
        self.files = files
        self.file_config = file_config
        self.op = op


class TB_operation_05B:
    """a class collect all config and function for temperature bias fit"""

    # initial guess and maxfev for the shared 2D temp-bias curve fit;
    # subclasses override when the shared default p0 does not converge
    TB_FIT_P0 = None
    TB_FIT_MAXFEV = 10000

    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        """config used in different function"""
        # universal config
        self.path = path
        self.files = [
            file
            for file in os.listdir(self.path)
            if "rundata" in file and "50C" not in file
        ]
        self.adc_max = 16384.0
        self.source = "Am241"
        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        file = os.path.join(self.path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="normal")
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def to_op(self):
        return Operation(self.files, lambda file: self.file_config(file), op=self)

    def load_data(self):
        data = (
            util.pickle_load(
                os.path.join(self.save_path, f"{os.path.splitext(file)[0]}.pickle")
            )
            for file in self.files
        )
        data_all = [[], [], [], []]
        for tb in data:
            fit_4ch = tb["fit_result"]
            tel_4ch = [
                {k: v[i] for k, v in tb["tel"].items() if len(v) == 4} for i in range(4)
            ]
            for data, fit, tel in zip(data_all, fit_4ch, tel_4ch):
                if fit is None:
                    continue
                center, center_err = fit["b"], fit["b_err"]
                temp, temp_err = np.average(tel["tempSipm"]), np.std(tel["tempSipm"])
                bias, bias_err = np.average(tel["bias"]), np.std(tel["bias"])
                data.append([center, center_err, temp, temp_err, bias, bias_err])
        data_all = [np.array(data) for data in data_all]
        return data_all

    def temp_bias_fit(self, data_all):
        result = []
        for ich, data in enumerate(data_all):
            center, center_err, temp, temp_err, bias, bias_err = (
                data[:, 0],
                data[:, 1],
                data[:, 2],
                data[:, 3],
                data[:, 4],
                data[:, 5],
            )
            try:
                res = util.temp_bias_fit_curvefit(
                    center,
                    center_err,
                    temp,
                    bias,
                    p0=self.TB_FIT_P0,
                    maxfev=self.TB_FIT_MAXFEV,
                )
            except util.FitError as e:
                print(f"chan {ich} fit failed: {e.args[-1]}")
                raise util.FitError(f"failed to do temp bias fit")
            result.append(res)

            data = np.stack([temp, bias], axis=1)
            name = f"temp_bias_fit_{ich}.png"
            plot.fit_err_plot_2d(
                data,
                center,
                lambda x: util.tempbias2DFunctionInternal(x, *(list(res.values())[:5])),
                ("temp$^\\circ$C", "bias/V", "center"),
                title=f"temp bias fit: channel {ich}",
                save_path=os.path.join(self.result_path, f"{util.headtime(name)}"),
            )
        util.json_save(
            result,
            os.path.join(self.result_path, f"{util.headtime('temp_bias_fit.json')}"),
        )
        return result


class EC_operation_05B:
    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        time_cut: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
        x_config,
    ) -> None:
        # read config
        self.x_config = x_config
        self.time_cut = util.json_load(time_cut)
        self.bkg_time_cut = {
            k: [v[1], v[2], v[3], v[0]] for k, v in self.time_cut.items()
        }

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 16384.0
        self.bin_width = 10
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_list = [
            f
            for f in os.listdir(self.x_path)
            if "_observe.dat" in f and "XM_22" not in f
        ]
        self.src_list = [
            f for f in os.listdir(self.src_path) if "rundata" in f and "bkg" not in f
        ]
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 52
        self.energy_split_low = 49

    def __get_src_bkg(self, file):
        src = file.split("_")[1]
        bkg = [
            f
            for f in os.listdir(self.src_path)
            if src in f and "bkg" in f and "rundata" in f
        ][0]
        return os.path.join(self.src_path, bkg)

    def xray_config(self, file):
        file = os.path.join(self.x_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(
            file,
            ending="xray",
            config_file=self.x_config,
            time_cut=self.time_cut[basename],
        )
        bkg_read_config = file_lib.Read_config(
            file,
            ending="xray",
            config_file=self.x_config,
            time_cut=self.bkg_time_cut[basename],
        )
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def src_config(self, file):
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        bkg = self.__get_src_bkg(basename)
        read_config = file_lib.Read_config(file, ending="normal")
        bkg_read_config = file_lib.Read_config(bkg, ending="normal")
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def to_x_op(self):
        return Operation(self.x_list, lambda file: self.xray_config(file), self)

    def to_src_op(self):
        return Operation(self.src_list, lambda file: self.src_config(file), self)

    def center_fit(self, energy: Float1D, center: Float1D, center_err: Float1D):
        popt, pcov = np.polyfit(
            center, energy, deg=2, full=False, cov=True, w=1.0 / center_err
        )
        perr = np.sqrt(np.diag(pcov))
        return list(popt), list(perr)

    def resolution_fit(self, energy, resolution, resolution_err):
        p0, pcov = util.resolution_polyfit(energy, resolution, resolution_err)
        perr = np.sqrt(np.diag(pcov))
        popt = list(p0)
        perr = list(perr)
        return popt, perr

    def ec_fit(self, src_result, src_energy, x_result, x_energy):
        result = src_result + x_result
        energy = src_energy + x_energy

        data = sorted(list(zip(energy, result)), key=lambda x: x[0])
        energy = np.array(list(map(lambda x: x[0], data)))
        center = [np.array([fit[i]["b"] for _, fit in data]) for i in range(4)]
        center_err = [np.array([fit[i]["b_err"] for _, fit in data]) for i in range(4)]
        resolution = [
            np.array([fit[i]["resolution"] for _, fit in data]) for i in range(4)
        ]
        resolution_err = [
            np.array([fit[i]["resolution_err"] for _, fit in data]) for i in range(4)
        ]
        q_low = energy < self.energy_split_low
        q_high = energy >= self.energy_split_high
        result = [{}, {}, {}, {}]
        for i in range(4):
            ec_low, ec_low_err = self.center_fit(
                energy[q_low], center[i][q_low], center_err[i][q_low]
            )
            ec_high, ec_high_err = self.center_fit(
                energy[q_high], center[i][q_high], center_err[i][q_high]
            )
            resolution_low, resolution_low_err = self.resolution_fit(
                energy[q_low], resolution[i][q_low], resolution_err[i][q_low]
            )
            resolution_high, resolution_high_err = self.resolution_fit(
                energy[q_high], resolution[i][q_high], resolution_err[i][q_high]
            )

            result[i] = {
                "channel": i,
                "EC_low": ec_low,
                "EC_low_err": ec_low_err,
                "EC_high": ec_high,
                "EC_high_err": ec_high_err,
                "resolution_low": resolution_low,
                "resolution_low_err": resolution_low_err,
                "resolution_high": resolution_high,
                "resolution_high_err": resolution_high_err,
            }
            fit_name = f"ec_coef_sci_ch{i}.json"
            util.json_save(result[i], f"{self.result_path}/{util.headtime(fit_name)}")
            save_data = np.array([energy, center[i]], dtype=np.float64)
            data_name = f"ec_data_ch{i}.npy"
            np.save(f"{self.result_path}/{util.headtime(data_name)}", arr=save_data)
        plot.ec_plot(
            energy,
            center,
            result,
            src_energy,
            x_energy,
            src_result,
            x_result,
            self.result_path,
            self.energy_split_low,
            self.energy_split_high,
        )
        return result


class TB_operation_03B(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        self.files = [
            file
            for file in os.listdir(self.path)
            if "rundata" in file
            and "baseline" not in file
            and "CI" not in file
            and "50C" not in file
        ]
        self.adc_max = 16384.0
        self.source = "Am241"

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        file = os.path.join(self.path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="03b")
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]


class EC_operation_03B(EC_operation_05B):
    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
    ) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 16384.0
        self.bin_width = 10
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_ch = [
            f for f in os.listdir(self.x_path) if "_rundata" in f and "CI" not in f
        ]

        self.x_list = list(set([f.split("_")[1] for f in self.x_ch]))
        self.x_list.sort()
        # self.src_list = [f for f in os.listdir(self.src_path) if 'rundata' in f and 'bkg' not in f and 'CI' not in f]
        self.src_list = [
            "src_Na22_20m_10cm_rundata2021-05-05-15-29-56.dat",
            "src_Am241_5m_10cm_rundata2021-05-05-15-15-48.dat",
            "src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat",
        ]
        self.src_bkg = [
            "src_bkg_5m_10cm_rundata2021-05-05-15-54-31.dat",
            "src_bkg_5m_10cm_rundata2021-05-05-14-55-58.dat",
            "src_bkg_5m_10cm_rundata2021-05-05-12-33-57.dat",
        ]
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 51
        self.energy_split_low = 49

    def resolution_fit(self, energy, resolution, resolution_err):
        popt, perr = util.resolution_lmfit(energy, resolution, resolution_err)
        return popt, perr

    def __get_x_files(self, energy_name: str):
        return [
            [
                os.path.join(self.x_path, f)
                for f in self.x_ch
                if f"{energy_name}_ch{i}" in f
            ][0]
            for i in range(4)
        ]

    def xray_config(self, energy_name: str):
        read_config = [
            file_lib.Read_config(ch_file, ending="03b", config_file=self.x_config)
            for ch_file in self.__get_x_files(energy_name)
        ]
        bkg_read_config = read_config[1:4]
        bkg_read_config.append(read_config[0])
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[energy_name], self.bkg_form[energy_name]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def __get_src_bkg(self, name: str):
        return self.src_bkg[self.src_list.index(name)]

    def src_config(self, file):
        bkg = os.path.join(self.src_path, self.__get_src_bkg(file))
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="03b-src")
        bkg_read_config = file_lib.Read_config(bkg, ending="03b-src")
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]


class TB_operation_07(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        #
        self.files = [
            file
            for file in os.listdir(self.path)
            if os.path.splitext(file)[1] == ".txt"
        ]
        self.adc_max = 65535.0
        self.source = "Am241"

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        file = os.path.join(self.path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="07")
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(
            bin_width=self.bin_width, adc_max=self.adc_max
        )
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]


class EC_operation_07(EC_operation_05B):
    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
    ) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 65535.0
        self.bin_width = 10
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_ch = [
            f
            for f in os.listdir(self.x_path)
            if "_ch" in f
            and "40p0" not in f
            and "15p0" not in f
            and "12p0" not in f
            and "99p9" not in f
            and "90p1" not in f
        ]

        self.x_list = list(set([f.split("_")[2] for f in self.x_ch]))
        self.x_list.sort()
        # self.src_list = [f for f in os.listdir(self.src_path) if 'rundata' in f and 'bkg' not in f and 'CI' not in f]
        self.src_list = [
            f
            for f in os.listdir(self.src_path)
            if "src" in f and "_bk_" not in f and "bkg" not in f
        ]
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 55
        self.energy_split_low = 49

    def __get_x_files(self, energy_name: str):
        return [
            [
                os.path.join(self.x_path, f)
                for f in self.x_ch
                if energy_name in f and f"_ch{i}_" in f
            ][0]
            for i in range(4)
        ]

    def xray_config(self, energy_name: str):
        read_config = [
            file_lib.Read_config(ch_file, ending="07", config_file=self.x_config)
            for ch_file in self.__get_x_files(energy_name)
        ]
        bkg_read_config = read_config[1:4]
        bkg_read_config.append(read_config[0])
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[energy_name], self.bkg_form[energy_name]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def __get_src_bkg(self, name: str):
        src = name.split("_")[3][:2]
        return [
            f for f in os.listdir(self.src_path) if "_bk_" not in f and f"bkg{src}" in f
        ][0]

    def src_config(self, file):
        bkg = os.path.join(self.src_path, self.__get_src_bkg(file))
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="07")
        bkg_read_config = file_lib.Read_config(bkg, ending="07")
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]


class TB_operation_04(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        self.files = [
            file
            for file in os.listdir(self.path)
            if os.path.splitext(file)[1] == ".txt"
        ]
        self.adc_max = 65535.0
        self.source = "Am241"

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        file = os.path.join(self.path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="04")
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(
            bin_width=self.bin_width, adc_max=self.adc_max
        )
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]


class EC_operation_04(EC_operation_05B):
    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
    ) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 65535.0
        self.bin_width = 10
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_ch = [
            f for f in os.listdir(self.x_path) if "_ch" in f and "_18p0_" not in f
        ]
        self.x_list = list(set([f.split("_")[3] for f in self.x_ch]))
        self.x_list.sort()

        self.src_list = [
            "210504162829_COM6_src_Co60_10m_10cm-Data.txt",
            "210504170604_COM6_src_Na22_30m_10cm-Data.txt",
            "210505120306_COM6_src_Cs137_15m_10cm-Data.txt",
            "210505151551_COM6_src_Am241_5m_10cm-Data.txt",
        ]
        self.src_bkg = [
            "210504164227_COM6_src_bkg_5m_10cm-Data.txt",
            "210504175217_COM6_src_bkg_5m_10cm-Data.txt",
            "210505123401_COM6_src_bkg_5m_10cm-Data.txt",
            "210505145547_COM6_src_bkg_5m_10cm-Data.txt",
        ]
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 55
        self.energy_split_low = 49

    def __get_x_files(self, energy_name: str):
        return [
            [
                os.path.join(self.x_path, f)
                for f in self.x_ch
                if f"{energy_name}_ch{i}" in f
            ][0]
            for i in range(4)
        ]

    def xray_config(self, energy_name: str):
        read_config = [
            file_lib.Read_config(ch_file, ending="04", config_file=self.x_config)
            for ch_file in self.__get_x_files(energy_name)
        ]
        bkg_read_config = read_config[1:4]
        bkg_read_config.append(read_config[0])
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[energy_name], self.bkg_form[energy_name]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def __get_src_bkg(self, name: str):
        return self.src_bkg[self.src_list.index(name)]

    def src_config(self, file):
        bkg = os.path.join(self.src_path, self.__get_src_bkg(file))
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="04")
        bkg_read_config = file_lib.Read_config(bkg, ending="04")
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def resolution_fit(self, energy, resolution, resolution_err):
        popt, perr = util.resolution_ExprFit(energy, resolution, resolution_err)
        return popt, perr


class TB_operation_10B(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        self.files = [
            file
            for file in os.listdir(self.path)
            if "observe" in file and "50C_265" not in file
        ]
        self.adc_max = 16384.0
        self.source = "Na22"

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        file = os.path.join(self.path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="10b")
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]


class EC_operation_10B(EC_operation_05B):
    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
    ) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 16384.0
        self.bin_width = 4
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_ch = [f for f in os.listdir(self.x_path) if "_ch" in f and "hk" not in f]
        self.x_list = list(set([f.split("_")[2] for f in self.x_ch]))
        self.x_list.sort()

        self.src_list = [
            "073_observe_Cs137.dat",
            "085_observe_Na22.dat",
            "089_observe_Am241.dat",
            "077_observe_Co60.dat",
        ]
        self.src_bkg = [
            "075_observe_Cs137_bkg.dat",
            "086_observe_Na22_bkg.dat",
            "090_observe_Am241_bkg.dat",
            "078_observe_Co60_bkg.dat",
        ]
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 55
        self.energy_split_low = 49

    def __get_x_files(self, energy_name: str):
        return [
            [
                os.path.join(self.x_path, f)
                for f in self.x_ch
                if f"{energy_name}_ch{i}" in f
            ][0]
            for i in range(4)
        ]

    def xray_config(self, energy_name: str):
        read_config = [
            file_lib.Read_config(ch_file, ending="10b")
            for ch_file in self.__get_x_files(energy_name)
        ]
        # bkg_read_config = read_config[1:4]
        # bkg_read_config.append(read_config[0])
        bkg_read_config = [
            read_config[1],
            read_config[2],
            read_config[0],
            read_config[0],
        ]
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[energy_name], self.bkg_form[energy_name]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def __get_src_bkg(self, name: str):
        return self.src_bkg[self.src_list.index(name)]

    def src_config(self, file):
        bkg = os.path.join(self.src_path, self.__get_src_bkg(file))
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="10b")
        bkg_read_config = file_lib.Read_config(bkg, ending="10b")
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def ec_fit(self, src_result, src_energy, x_result, x_energy):
        result = src_result + x_result
        energy = src_energy + x_energy
        CHN_NUM = 3
        data = sorted(list(zip(energy, result)), key=lambda x: x[0])
        energy = np.array(list(map(lambda x: x[0], data)))
        center = [np.array([fit[i]["b"] for _, fit in data]) for i in range(CHN_NUM)]
        center_err = [
            np.array([fit[i]["b_err"] for _, fit in data]) for i in range(CHN_NUM)
        ]
        resolution = [
            np.array([fit[i]["resolution"] for _, fit in data]) for i in range(CHN_NUM)
        ]
        resolution_err = [
            np.array([fit[i]["resolution_err"] for _, fit in data])
            for i in range(CHN_NUM)
        ]
        q_low = energy < self.energy_split_low
        q_high = energy >= self.energy_split_high
        result = [{}, {}, {}, {}]
        for i in range(CHN_NUM):
            ec_low, ec_low_err = self.center_fit(
                energy[q_low], center[i][q_low], center_err[i][q_low]
            )
            ec_high, ec_high_err = self.center_fit(
                energy[q_high], center[i][q_high], center_err[i][q_high]
            )
            resolution_low, resolution_low_err = self.resolution_fit(
                energy[q_low], resolution[i][q_low], resolution_err[i][q_low]
            )
            resolution_high, resolution_high_err = self.resolution_fit(
                energy[q_high], resolution[i][q_high], resolution_err[i][q_high]
            )

            result[i] = {
                "channel": i,
                "EC_low": ec_low,
                "EC_low_err": ec_low_err,
                "EC_high": ec_high,
                "EC_high_err": ec_high_err,
                "resolution_low": resolution_low,
                "resolution_low_err": resolution_low_err,
                "resolution_high": resolution_high,
                "resolution_high_err": resolution_high_err,
            }
            fit_name = f"ec_coef_sci_ch{i}.json"
            util.json_save(result[i], f"{self.result_path}/{util.headtime(fit_name)}")
            save_data = np.array([energy, center[i]], dtype=np.float64)
            data_name = f"ec_data_ch{i}.npy"
            np.save(f"{self.result_path}/{util.headtime(data_name)}", arr=save_data)
        # 临时补丁，让数据为四通道
        center = [center[0], center[1], center[2], center[0]]
        result = [result[0], result[1], result[2], result[0]]
        src_result = [[s[0], s[1], s[2], s[0]] for s in src_result]
        x_result = [[s[0], s[1], s[2], s[0]] for s in x_result]
        plot.ec_plot(
            energy,
            center,
            result,
            src_energy,
            x_energy,
            src_result,
            x_result,
            self.result_path,
            self.energy_split_low,
            self.energy_split_high,
        )
        return result


class TB_operation_11B(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        # 全路径
        self.tb_files = self.get_tb_files(self.path)
        # 文件名
        self.files = [Path(f).name for f in self.tb_files]
        self.adc_max = 16384.0
        self.source = "Cs137"

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def get_tb_files(self, path):
        path = Path(path)
        paths = [
            path / r"温度-偏压实验-20~-10℃",
            path / r"温度偏压0～10摄氏度",
            path / r"温度偏压20～30摄氏度",
            path / r"温度偏压40～50摄氏度",
        ]
        tb_files = [list(Path(p).glob("*_observe*.dat")) for p in paths]
        tb_files = [
            str(item)
            for sublist in tb_files
            for item in sublist
            if util.not_contain(item, "on", "off")
        ]
        tb_files.remove(
            str(path / r"温度偏压0～10摄氏度/247_0_Cs137_27.5_observe_1.dat")
        )
        tb_files.remove(
            str(path / r"温度偏压0～10摄氏度/004_10_Cs137_26.5_1_observe.dat")
        )
        tb_files.remove(str(path / r"温度偏压40～50摄氏度/039_Cs_40_26.5_observe.dat"))
        tb_files.remove(str(path / r"温度-偏压实验-20~-10℃/232_Cs_-10_observe.dat"))
        tb_files.remove(
            str(path / r"温度-偏压实验-20~-10℃/233_Cs_-10_26.5_observe.dat")
        )
        tb_files = [f for f in tb_files if "_50_Cs_2" not in f]
        tb_files.remove(str(path / r"温度偏压0～10摄氏度/013_10_Cs137_29_observe .dat"))
        tb_files.sort()

        return tb_files

    def get_key(self, file: str):
        return Path(file).stem

    def get_path(self, file: str):
        key = self.get_key(file)
        files = [f for f in self.tb_files if self.get_key(f) == key]
        if len(files) == 0:
            raise FileNotFoundError(f"File with key {key} not found.")
        return files[0]

    def file_config(self, file):
        path = self.get_path(file)
        key = self.get_key(file)
        read_config = file_lib.Read_config(path, ending="11b")
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(
            bin_width=self.bin_width, adc_max=self.adc_max
        )
        fit_config = file_lib.Fit_config(self.fit_range[key])
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def temp_bias_fit(self, data_all):
        result = []
        for ich, data in enumerate(data_all):
            center, center_err, temp, temp_err, bias, bias_err = (
                data[:, 0],
                data[:, 1],
                data[:, 2],
                data[:, 3],
                data[:, 4],
                data[:, 5],
            )
            try:
                res = util.temp_bias_lmfit(
                    center, center_err, temp, temp_err, bias, bias_err
                )
            except util.FitError as e:
                print(f"chan {ich} fit failed: {e.args[-1]}")
                raise util.FitError(f"failed to do temp bias fit")
            result.append(res)

            data = np.stack([temp, bias], axis=1)
            name = f"temp_bias_fit_{ich}.png"
            plot.fit_err_plot_2d(
                data,
                center,
                lambda x: util.tempbias2DFunctionInternal(x, *(list(res.values())[:5])),
                ("temp$^\\circ$C", "bias/V", "center"),
                title=f"temp bias fit: channel {ich}",
                save_path=os.path.join(self.result_path, f"{util.headtime(name)}"),
            )
        util.json_save(
            result,
            os.path.join(self.result_path, f"{util.headtime('temp_bias_fit.json')}"),
        )
        return result


class EC_operation_11B(EC_operation_05B):
    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
    ) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 16384.0
        self.bin_width = 4
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_ch = [f for f in os.listdir(self.x_path) if "_ch" in f and "hk" not in f]
        self.x_list = list(set([f.split("_")[2] for f in self.x_ch]))
        # 道址异常高，暂时去除
        self.x_list.remove("20")
        self.x_list.sort()

        self.src_list = [
            "182_Ba133_20min_observe.dat",
            "187_Cs137_ch012_3min_observe.dat",
            "189_Eu152_2min_observe.dat",
            "191_Co60_2min_observe.dat",
        ]
        self.src_bkg = [
            "",
            "188_Cs137_ch3_15min_observe.dat",
            "190_Eu152_15min_ch3_observe.dat",
            "192_Co60_15min_ch3_observe.dat",
        ]
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 55
        self.energy_split_low = 49

    def __get_x_files(self, energy_name: str):
        return [
            [
                os.path.join(self.x_path, f)
                for f in self.x_ch
                if f"{energy_name}_ch{i}" in f
            ][0]
            for i in range(4)
        ]

    def xray_config(self, energy_name: str):
        read_config = [
            file_lib.Read_config(ch_file, ending="11b")
            for ch_file in self.__get_x_files(energy_name)
        ]
        # bkg_read_config = read_config[1:4]
        # bkg_read_config.append(read_config[0])
        bkg_read_config = [
            read_config[1],
            read_config[2],
            read_config[0],
            read_config[0],
        ]
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[energy_name], self.bkg_form[energy_name]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def __get_src_bkg(self, name: str):
        return self.src_bkg[self.src_list.index(name)]

    def src_config(self, file):
        bkg_name = self.__get_src_bkg(file)
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="11b")
        if bkg_name == "":
            bkg_read_config = file_lib.Read_config("", ending="11b")
        else:
            bkg_read_config = file_lib.Read_config(
                os.path.join(self.src_path, bkg_name), ending="11b"
            )
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def ec_fit(self, src_result, src_energy, x_result, x_energy):
        result = src_result + x_result
        energy = src_energy + x_energy
        CHN_NUM = 3
        data = sorted(list(zip(energy, result)), key=lambda x: x[0])
        energy = np.array(list(map(lambda x: x[0], data)))
        center = [np.array([fit[i]["b"] for _, fit in data]) for i in range(CHN_NUM)]
        center_err = [
            np.array([fit[i]["b_err"] for _, fit in data]) for i in range(CHN_NUM)
        ]
        resolution = [
            np.array([fit[i]["resolution"] for _, fit in data]) for i in range(CHN_NUM)
        ]
        resolution_err = [
            np.array([fit[i]["resolution_err"] for _, fit in data])
            for i in range(CHN_NUM)
        ]
        q_low = energy < self.energy_split_low
        q_high = energy >= self.energy_split_high
        result = [{}, {}, {}, {}]
        for i in range(CHN_NUM):
            ec_low, ec_low_err = self.center_fit(
                energy[q_low], center[i][q_low], center_err[i][q_low]
            )
            ec_high, ec_high_err = self.center_fit(
                energy[q_high], center[i][q_high], center_err[i][q_high]
            )
            resolution_low, resolution_low_err = self.resolution_fit(
                energy[q_low], resolution[i][q_low], resolution_err[i][q_low]
            )
            resolution_high, resolution_high_err = self.resolution_fit(
                energy[q_high], resolution[i][q_high], resolution_err[i][q_high]
            )

            result[i] = {
                "channel": i,
                "EC_low": ec_low,
                "EC_low_err": ec_low_err,
                "EC_high": ec_high,
                "EC_high_err": ec_high_err,
                "resolution_low": resolution_low,
                "resolution_low_err": resolution_low_err,
                "resolution_high": resolution_high,
                "resolution_high_err": resolution_high_err,
            }
            fit_name = f"ec_coef_sci_ch{i}.json"
            util.json_save(result[i], f"{self.result_path}/{util.headtime(fit_name)}")
            save_data = np.array([energy, center[i]], dtype=np.float64)
            data_name = f"ec_data_ch{i}.npy"
            np.save(f"{self.result_path}/{util.headtime(data_name)}", arr=save_data)
        # 临时补丁，让数据为四通道
        center = [center[0], center[1], center[2], center[0]]
        result = [result[0], result[1], result[2], result[0]]
        src_result = [[s[0], s[1], s[2], s[0]] for s in src_result]
        x_result = [[s[0], s[1], s[2], s[0]] for s in x_result]
        plot.ec_plot(
            energy,
            center,
            result,
            src_energy,
            x_energy,
            src_result,
            x_result,
            self.result_path,
            self.energy_split_low,
            self.energy_split_high,
        )
        return result

class TB_operation_09(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        self.files = [
            file
            for file in os.listdir(self.path)
            if os.path.splitext(file)[1] == ".txt"
        ]
        self.files.remove("0826_30C_265_5m_0x0090.txt")
        self.files.remove("0826_30C_265_5m_0x00A0.txt")
        self.files.remove("0826_30C_265_5m_0x00B0.txt")
        self.files.remove("0826_30C_265_5m_0x00BB.txt")
        self.files.remove("0826_20C_290_5m_0x01F0 (1).txt")
        self.files.remove("0825_0C_290_4m_0x00C5.txt")

        
        self.adc_max = 65535.0
        self.source = "Am241"

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        file = os.path.join(self.path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="09")
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(
            bin_width=self.bin_width, adc_max=self.adc_max
        )
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]
    
class EC_operation_09(EC_operation_05B):
    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
    ) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 65535.0
        self.bin_width = 10
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_ch = [
            f for f in os.listdir(self.x_path) if "_ch" in f and "65keV_" not in f
        ]
        self.x_list = list(set([f.split("_")[0] for f in self.x_ch]))
        self.x_list.sort()

        self.src_list = [
            "0826_10C_285_Na22_20m_0x00CF.txt",
            "0827_10C_285_Co60_20m_0x010F.txt",
            "0827_10C_285_Cs137_20m_0x010F.txt"

        ]
        self.src_bkg = [
            "0826_10C_285_Na22_bkg_20m.txt",
            "0827_10C_285_bkg_20m_0x00CF.txt",
            "0827_10C_285_bkg_20m_0x00CF.txt",
        ]
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 55
        self.energy_split_low = 49

    def __get_x_files(self, energy_name: str):
        return [
            [
                os.path.join(self.x_path, f)
                for f in self.x_ch
                if f"{energy_name}_ch{i}" in f
            ][0]
            for i in range(4)
        ]

    def xray_config(self, energy_name: str):
        read_config = [
            file_lib.Read_config(ch_file, ending="09", config_file=self.x_config)
            for ch_file in self.__get_x_files(energy_name)
        ]
        bkg_read_config = read_config[1:4]
        bkg_read_config.append(read_config[0])
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range.get(energy_name, [[None, None]]*4), self.bkg_form.get(energy_name, "lin")
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def __get_src_bkg(self, name: str):
        return self.src_bkg[self.src_list.index(name)]

    def src_config(self, file):
        bkg = os.path.join(self.src_path, self.__get_src_bkg(file))
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="09")
        bkg_read_config = file_lib.Read_config(bkg, ending="09")
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range.get(basename, [[None, None]]*4), self.bkg_form.get(basename, "lin")
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

class TB_operation_12B(TB_operation_05B):
    """12B TB data read from the backup directory per tb_file_map.json
    (see docs/12B_13B/data.md): all 7 temperature x 8 bias points, each
    mapped to a 备份/072/sample/{NNN}_observe.dat + hk/ecu_1_{NNN}.hk pair.

    Processing names are "{temp}C_{bias}" (e.g. "-10C_265", "30C_290").
    """

    # points excluded (with reason):
    # - (-20, 275): hk ecu_1_060 has two bias epochs (26.1 V then 28.0 V) and
    #   the sci data belongs to the 26.1 V epoch while the plateau is 28.0 V
    # - (-20, 285): temperature drifts -11.3 -> -7.7 C within the run
    #   (Tspread 0.9 C vs ~0.03 C for normal files) and bias sags at the end
    EXCLUDE = [(-20, 275), (-20, 285)]

    def __init__(self, path, fit_range, save_path, save_fig_path, result_path, file_map) -> None:
        self.path = path  # base dir for the map's relative paths (raw_data root)
        points = util.json_load(file_map)
        self.point_map = {}
        for p in points:
            if (p["temp_setpoint_C"], p["bias_code"]) in self.EXCLUDE:
                continue
            name = f"{p['temp_setpoint_C']}C_{p['bias_code']}"
            self.point_map[name] = p
        self.files = sorted(self.point_map.keys())
        self.adc_max = 16384.0
        self.source = "Am241"

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        p = self.point_map[file]
        sci_path = os.path.join(self.path, p["observe_file"])
        hk_path = os.path.join(self.path, p["hk_file"])
        kwarg = {"hk_path": hk_path}
        # shared observe file: -10C/270 takes the first half (27.0 V epoch),
        # -10C/275 the second half (27.5 V epoch)
        if "前半段" in p["note"]:
            kwarg.update(sci_half="first", hk_bias=p["bias_setpoint_V"])
        elif "后半段" in p["note"]:
            kwarg.update(sci_half="second", hk_bias=p["bias_setpoint_V"])
        read_config = file_lib.Read_config(sci_path, ending="12b", kwarg=kwarg)
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(
            bin_width=self.bin_width, adc_max=self.adc_max
        )
        fit_config = file_lib.Fit_config(self.fit_range[file])
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    # 12B-appropriate initial guess: the shared default p0 does not converge
    # for 12B, so the base-class temp_bias_fit uses these instead
    TB_FIT_P0 = [-0.02, 0.07, 24.4, -35.0, -1000.0]
    TB_FIT_MAXFEV = 100000

    def load_data(self):
        data_all = super().load_data()
        # the Vov^2 response model under-predicts the gain at the lowest
        # fitted bias (27.0 V row, residuals +3~+6% while all other points
        # are < 3%): the shared model form cannot follow the response there,
        # so the 2D fit is restricted to bias >= 27.5 V. The 27.0 V single
        # fits are still produced and kept in the pickles.
        return [data[data[:, 4] >= 27.25] for data in data_all]


class EC_operation_12B(EC_operation_05B):
    """12B EC data: per-channel X-ray files {idx}_{kv}_ch{n}.dat grouped by kV,
    4 source files with a shared environmental background 0611env.dat."""

    def __init__(
        self,
        tb_result_path: str,
        fit_range: str,
        energy: str,
        bkg_form: str,
        x_path,
        src_path,
        save_path,
        save_fig_path,
        result_path,
    ) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5
        tb_result: List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [
            lambda t, b: util.tempbias2DFunction(
                t, b, c["G0"], c["k"], c["V0"], c["b"], c["c"]
            )
            for c in tb_result
        ]
        self.corr = [lambda t, b: f(ref_temp, ref_bias) / f(t, b) for f in ref_func]
        self.adc_max = 16384.0
        self.bin_width = 4
        # fit config
        self.fit_range = util.json_load(fit_range)
        self.bkg_form = util.json_load(bkg_form)
        # ec file process
        self.x_path = x_path
        self.src_path = src_path
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
        self.x_ch = [
            f
            for f in os.listdir(self.x_path)
            if f.endswith(".dat") and "_ch" in f and "old" not in f
        ]
        self.x_list = list(set([f.split("_")[1] for f in self.x_ch]))
        self.x_list.sort(key=int)
        # drop energy points whose hk pairing is incomplete (e.g. 90 kV has
        # only _old hk, 75 kV lacks ch1 hk)
        self.x_list = [e for e in self.x_list if self.__x_hk_complete(e)]
        # drop energy points without a full set of 4-channel fit ranges:
        # the shared split ec_fit / plot.ec_plot require every point to have
        # all 4 channels (points with a null channel, e.g. 15 kV ch3, must
        # be dropped as a whole)
        self.x_list = [
            e
            for e in self.x_list
            if e in self.fit_range and all(r is not None for r in self.fit_range[e])
        ]

        # Am241 run: hk mean bias 20.4 V is a ramp-artifact -- the hk has a
        # 422/604-record 28.5 V plateau (the dominant epoch) and the sci
        # spectrum shows the 59.5 keV line, so it is a valid EC point.
        # 0611env.dat is the shared environmental background.
        self.src_list = [
            "0611_Am241_12min_240f0032.dat",
            "0611_Na22_10min_240f0032.dat",
            "0611_Cs137_30min_240f0064.dat",
            "0611_Co60_25min_240f0032.dat",
        ]
        self.src_bkg = ["0611env.dat"] * 4
        self.energy = util.json_load(energy)
        # self.energy_split = 50.2 # keV, absorption edges of Gd
        self.energy_split_high = 55
        self.energy_split_low = 49

    def __get_x_files(self, energy_name: str):
        return [
            [
                os.path.join(self.x_path, f)
                for f in self.x_ch
                if f.split("_")[1] == energy_name and f"_ch{i}" in f
            ][0]
            for i in range(4)
        ]

    def __x_hk_complete(self, energy_name: str) -> bool:
        from lib_reader.reader12.read import getHK

        try:
            files = self.__get_x_files(energy_name)
            for f in files:
                getHK(f)
        except (IndexError, FileNotFoundError):
            return False
        return True

    def xray_config(self, energy_name: str):
        read_config = [
            file_lib.Read_config(ch_file, ending="12b")
            for ch_file in self.__get_x_files(energy_name)
        ]
        bkg_read_config = [
            read_config[1],
            read_config[2],
            read_config[0],
            read_config[0],
        ]
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[energy_name], self.bkg_form[energy_name]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    def __get_src_bkg(self, name: str):
        return self.src_bkg[self.src_list.index(name)]

    def src_config(self, file):
        bkg = os.path.join(self.src_path, self.__get_src_bkg(file))
        file = os.path.join(self.src_path, file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending="12b")
        bkg_read_config = file_lib.Read_config(bkg, ending="12b")
        spectrum_config = file_lib.Spectrum_config(
            corr=self.corr, adc_max=self.adc_max, bin_width=self.bin_width
        )
        fit_config = file_lib.Fit_config(
            self.fit_range[basename], self.bkg_form[basename]
        )
        return [read_config, bkg_read_config, spectrum_config, fit_config]


class TB_operation_N1(TB_operation_05B):
    """N1 (GRIDN1) TB data: GAGG runs (584-byte ft packets) and CLYC runs
    (1080-byte wf packets) under {path}/GAGG and {path}/CLYC.

    GAGG negative temperatures are single-bias files; all other runs (and all
    CLYC runs) are bias-scan files split into per-bias segments by the reader
    (seg_bias kwarg). Files with "To" (temperature transition), "CI" (current
    scan) or "test" in the name are excluded, mirroring the dataset selection
    of the gridN_cali L0 notebooks.

    Processing names are "{ds}_{temp}_{bias}" (e.g. "GAGG_m20C_265").

    Sources differ per dataset: GAGG runs used an Am241 source and cover
    channels 1/2; CLYC runs used a Na22 source and cover all 4 channels.
    """

    # per-dataset source annotation (the two datasets use different sources)
    source = {"GAGG": "Am241", "CLYC": "Na22"}

    # CLYC scans carry extra re-measured bias segments (same per-file list as
    # the gridN_cali read_raw_CLYC.ipynb vol_sets)
    CLYC_EXTRA_BIAS = {
        "m20C": (277, 279, 281),
        "m10C": (279, 281),
        "0C": (281,),
        "10C": (289,),
        "20C": (289, 291),
        "30C": (289, 291),
    }

    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        import re

        self.path = path
        self.point_map = {}
        for ds, mode in (("GAGG", "ft"), ("CLYC", "wf")):
            d = os.path.join(path, ds)
            for f in sorted(os.listdir(d)):
                if not f.endswith(".event.dat") or "To" in f or "CI" in f or "test" in f:
                    continue
                m = re.match(r"^(m?\d+C)-(\d{3})-\d+\.event\.dat$", f)
                if m:  # single-bias point
                    temp, bias = m.group(1), int(m.group(2))
                    self.point_map[f"{ds}_{temp}_{bias}"] = (os.path.join(d, f), mode, None)
                    continue
                m = re.match(r"^(m?\d+C)-265-290-\d+\.event\.dat$", f)
                if m:  # bias scan: one segment per bias code
                    temp = m.group(1)
                    biases = [265, 270, 275, 280, 283, 285, 287, 290]
                    if ds == "CLYC":
                        biases += list(self.CLYC_EXTRA_BIAS.get(temp, ()))
                    for bias in biases:
                        self.point_map[f"{ds}_{temp}_{bias}"] = (
                            os.path.join(d, f), mode, bias,
                        )
        self.files = sorted(self.point_map.keys())
        self.adc_max = 16384.0

        self.fit_range = util.json_load(fit_range)
        self.bin_width = 4
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path

    def file_config(self, file):
        sci_path, mode, seg_bias = self.point_map[file]
        kwarg = {"mode": mode}
        if seg_bias is not None:
            kwarg["seg_bias"] = seg_bias
        read_config = file_lib.Read_config(sci_path, ending="n1", kwarg=kwarg)
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(
            bin_width=self.bin_width, adc_max=self.adc_max
        )
        fit_config = file_lib.Fit_config(self.fit_range[file])
        return [read_config, bkg_read_config, spectrum_config, fit_config]

    # per-dataset, per-channel initial guesses for the 2D fit, taken from the
    # gridN_cali L2 results (same 5-parameter model); the two datasets were
    # taken with different sources (GAGG: Am241, CLYC: Na22) so their peak
    # positions are on different ADC scales and are fitted separately
    TB_FIT_P0_BY_DS = {
        "GAGG": {
            1: [3.296e-04, 0.0182, 23.87, 64.7, 5.197e4],
            2: [3.775e-04, 0.0188, 23.98, 66.6, 4.476e4],
        },
        "CLYC": {
            0: [7.305e-04, 0.0162, 23.60, 500.0, 3.201e4],
            1: [1.062e-02, 0.0195, 23.81, 49.7, 1.561e4],
            2: [1.072e-02, 0.0191, 23.97, 47.1, 1.560e4],
            3: [7.776e-04, 0.0129, 23.57, 500.0, 3.721e4],
        },
    }

    # which dataset each channel is fitted on (GAGG runs target ch1/2,
    # CLYC runs cover all 4 channels)
    CHANNEL_DS = {0: "CLYC", 1: "GAGG", 2: "GAGG", 3: "CLYC"}

    def load_data(self):
        """like the base class but grouped per dataset (GAGG/CLYC, the two DAQ
        modes are on different ADC scales) and skipping qa_flag == 'fail'
        single fits (redchi above fail threshold)"""
        data_by_ds = {"GAGG": [[], [], [], []], "CLYC": [[], [], [], []]}
        for file in self.files:
            tb = util.pickle_load(
                os.path.join(self.save_path, f"{os.path.splitext(file)[0]}.pickle")
            )
            ds = "GAGG" if "/GAGG/" in tb["file"] else "CLYC"
            fit_4ch = tb["fit_result"]
            tel_4ch = [
                {k: v[i] for k, v in tb["tel"].items() if len(v) == 4} for i in range(4)
            ]
            for data, fit, tel in zip(data_by_ds[ds], fit_4ch, tel_4ch):
                if fit is None or fit.get("qa_flag") == "fail":
                    continue
                center, center_err = fit["b"], fit["b_err"]
                temp, temp_err = np.average(tel["tempSipm"]), np.std(tel["tempSipm"])
                bias, bias_err = np.average(tel["bias"]), np.std(tel["bias"])
                data.append([center, center_err, temp, temp_err, bias, bias_err])
        return {ds: [np.array(d) for d in v] for ds, v in data_by_ds.items()}

    def temp_bias_fit(self, data_by_ds):
        """fit each dataset separately (the two sources put the peaks on
        different ADC scales), write the two parts
        (temp_bias_fit_am241.json = GAGG ch1/2, temp_bias_fit_na22_511.json =
        CLYC all 4 channels; both 4-slot arrays with nulls) and the standard
        merged temp_bias_fit.json (ch1/2 from GAGG, ch0/3 from CLYC)."""
        part_ds = {"am241": "GAGG", "na22_511": "CLYC"}
        part_result = {}
        for part, ds in part_ds.items():
            data_all = data_by_ds[ds]
            result = [None] * 4
            for ich, data in enumerate(data_all):
                if ich not in self.TB_FIT_P0_BY_DS[ds] or len(data) < 6:
                    continue
                center, center_err, temp, bias = (
                    data[:, 0], data[:, 1], data[:, 2], data[:, 4],
                )
                try:
                    res = util.temp_bias_fit_curvefit(
                        center, center_err, temp, bias,
                        p0=self.TB_FIT_P0_BY_DS[ds][ich], maxfev=100000,
                    )
                except util.FitError as e:
                    print(f"{ds} chan {ich} fit failed: {e.args[-1]}")
                    continue
                result[ich] = res
                # canonical per-channel plots come from the merged source
                if self.CHANNEL_DS[ich] != ds:
                    continue
                plot.fit_err_plot_2d(
                    np.stack([temp, bias], axis=1),
                    center,
                    lambda x: util.tempbias2DFunctionInternal(x, *(list(res.values())[:5])),
                    ("temp$^\\circ$C", "bias/V", "center"),
                    title=f"temp bias fit: channel {ich} ({ds})",
                    save_path=os.path.join(
                        self.result_path, f"{util.headtime('temp_bias_fit_' + str(ich) + '.png')}"
                    ),
                )
            util.json_save(
                result,
                os.path.join(
                    self.result_path, f"{util.headtime(f'temp_bias_fit_{part}.json')}"
                ),
            )
            part_result[part] = result

        merged = [None] * 4
        for ich, ds in self.CHANNEL_DS.items():
            part = "am241" if ds == "GAGG" else "na22_511"
            merged[ich] = part_result[part][ich]
        util.json_save(
            merged,
            os.path.join(self.result_path, f"{util.headtime('temp_bias_fit.json')}"),
        )
        return merged


class EC_operation_N1(EC_operation_05B):
    """N1 EC: not onboarded yet (GRIDN1 EC data live in the sibling
    experiment directories 260326/260327放射源 and 260129/260202计量院标定).
    This stub only stores the config so that `just tb N1 ...` works; the EC
    processing will be implemented when the EC data are organized."""

    def __init__(self, **kwargs) -> None:
        self.__dict__.update(kwargs)

    def to_x_op(self):
        return Operation([], lambda file: None, self)

    def to_src_op(self):
        return Operation([], lambda file: None, self)


def __get_fp05B(config, nocache=False) -> file_lib.File_operation_05b:
    return file_lib.File_operation_05b(config[0].path, *config, nocache=nocache)


def __dict_4ch_reconstruct(dict_4ch):
    # dict_4ch is list of 4 dict all value to be list with 4 elements, return a dict with all value to be 4 ch, each value from corresponding dict
    dict_4ch_re = {}
    for key in dict_4ch[0].keys():
        dict_4ch_re[key] = []
        for i in range(4):
            if len(dict_4ch[i][key]) == 4:
                dict_4ch_re[key].append(dict_4ch[i][key][i])
            else:
                dict_4ch_re[key].append(dict_4ch[i][key])
    return dict_4ch_re


def __get_fp03B(config, nocache=False) -> file_lib.File_operation_05b:
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fps = [
        file_lib.File_operation_05b(
            read_config[i].path,
            read_config[i],
            bkg_read_config[i],
            spectrum_config,
            fit_config,
            nocache=nocache,
        )
        for i in range(4)
    ]
    sci = __dict_4ch_reconstruct([fps[i].sci for i in range(4)])
    tel = __dict_4ch_reconstruct([fps[i].tel for i in range(4)])
    bkg_sci = __dict_4ch_reconstruct([fps[i].bkg_sci for i in range(4)])
    bkg_tel = __dict_4ch_reconstruct([fps[i].bkg_tel for i in range(4)])
    fp = fps[0]
    fp.sci, fp.tel = sci, tel
    fp.bkg_sci, fp.bkg_tel = bkg_sci, bkg_tel
    return fp


def process(op: Operation, file: str, fp_method=None, **kw_args) -> None:
    """process file in sinlge experiment, like specific temp bias or photon energy

    Parameters
    ----------
    op: Operation
        from TB_operation.to_op, EC_operation.to_x_op, EC_operation.to_src_op
    file: str
        element of op.files
    fp_method: Optional, default to be None for 05B like, else to be 03B like

    kw_args:
        x_lim: [int, int], optional, default to be None
            plot the raw figure, specify the x_lim
        save_path: str, optional, default to be None
            path of the figure to save, None for show the figure
    """
    fp_method = __get_fp05B if fp_method is None else __get_fp03B
    config = op.file_config(file)
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fp = fp_method(config, nocache=kw_args.get("nocache", False))
    # inject QA thresholds so peak_fit can classify results
    save_path = op.op.save_path
    ver = Path(save_path).parts[1]  # e.g. "data/03B/single_process/..." -> "03B"
    category = "tb" if "TB_fit" in save_path else "ec"
    fp.qa_thresholds = util.load_qa_thresholds(ver, category)
    fp.get_spectrum()
    # plot raw spectrum
    if kw_args.get("x_lim", None) is not None:
        plot.raw_plot(
            fp.spectrum,
            fp.x,
            title=file,
            x_lim=kw_args["x_lim"],
            save_path=kw_args.get("save_path", None),
        )
        return
    # fit and save fit result
    fp.peak_fit()
    file = os.path.splitext(file)[0]
    plot.fit_plot(
        fp.spectrum,
        fp.x,
        fp.fit_result,
        title=f"{file}: {np.mean(fp.tel['bias'][0]):.2f}V, {np.mean(fp.tel['tempSipm'][0]):.2f}C",
        bkgForm=fit_config.bkg_form,
        fit_range=fit_config.fit_range,
        save_path=f"{op.op.save_fig_path}/{file}.png",
    )
    fp.save(os.path.join(op.op.save_path, f"{file}.pickle"))
