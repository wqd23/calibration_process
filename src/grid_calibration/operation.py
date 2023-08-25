import lib_reader as ver
from . import util_lib as util
from lib_reader.reader05.my_type import *
import os
import numpy as np
from dataclasses import dataclass
import json
import matplotlib.pyplot as plt
from . import file_lib


class Operation():
    '''wrapper for tb_op, x of ec_op, and src of ec_op
    '''
    def __init__(self, files, file_config, op) -> None:
        self.files = files
        self.file_config = file_config
        self.op = op
class TB_operation_05B():
    """a class collect all config and function for temperature bias fit
    """
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        """config used in different function
        """
        # universal config
        self.path = path
        self.files = [file for file in os.listdir(self.path) if "rundata" in file and '50C' not in file]
        self.adc_max = 16384.0
        self.source = "Am241"
        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
    def file_config(self, file):
        file = os.path.join(self.path,file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending='normal')
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]
    def to_op(self):
        return Operation(self.files, lambda file: self.file_config(file), op=self)
    def load_data(self):
        data = (util.pickle_load(os.path.join(self.save_path, f"{os.path.splitext(file)[0]}.pickle")) for file in self.files)
        data_all = [[],[],[],[]]
        for tb in data:
            fit_4ch = tb['fit_result']
            tel_4ch = [{k:v[i] for k,v in tb['tel'].items() if len(v) == 4} for i in range(4)]
            for data, fit, tel in zip(data_all, fit_4ch, tel_4ch):
                if fit is None:
                    continue
                center, center_err = fit['b'], fit['b_err']
                temp, temp_err = np.average(tel['tempSipm']), np.std(tel['tempSipm'])
                bias, bias_err = np.average(tel['bias']), np.std(tel['bias'])
                data.append([center, center_err, temp, temp_err, bias, bias_err])
        data_all = [np.array(data) for data in data_all]
        return data_all

    def temp_bias_fit(self,data_all):
        result = []
        for ich,data in enumerate(data_all):
            center, center_err, temp, temp_err, bias, bias_err = data[:,0],data[:,1],data[:,2],data[:,3],data[:,4],data[:,5]
            try:
                res = util.temp_bias_fit_curvefit(center, center_err, temp, bias)
            except util.FitError as e:
                print(f"chan {ich} fit failed: {e.args[-1]}")
                raise util.FitError(f"failed to do temp bias fit")
            result.append(res)

            data = np.stack([temp,bias], axis=1)
            name = f'temp_bias_fit_{ich}.png'
            util.fit_err_plot_2d(data, center, lambda x: util.tempbias2DFunctionInternal(x, *(list(res.values())[:5])),("temp$^\\circ$C","bias/V","center"), title=f"temp bias fit: channel {ich}", save_path=os.path.join(self.result_path, f"{util.headtime(name)}"))
        util.json_save(result, os.path.join(self.result_path, f"{util.headtime('temp_bias_fit.json')}"))
        return result

class EC_operation_05B():
    def __init__(self, tb_result_path:str, fit_range:str, time_cut:str, energy:str, bkg_form:str,x_path, src_path, save_path, save_fig_path, result_path, x_config) -> None:
        # read config
        self.x_config = x_config
        self.time_cut = util.json_load(time_cut)
        self.bkg_time_cut = {k:[v[1],v[2],v[3],v[0]] for k,v in self.time_cut.items()}

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5 
        tb_result:List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [lambda t,b: util.tempbias2DFunction(t, b, c['G0'],c['k'],c['V0'],c['b'],c['c']) for c in tb_result]
        self.corr = [lambda t,b: f(ref_temp, ref_bias)/f(t,b) for f in ref_func]
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
        self.x_list = [f for f in os.listdir(self.x_path) if '_observe.dat' in f and 'XM_22' not in f]
        self.src_list = [f for f in os.listdir(self.src_path) if 'rundata' in f and 'bkg' not in f]
        self.energy = util.json_load(energy)
        self.energy_split = 50.2 # keV, absorption edges of Gd
    def __get_src_bkg(self, file):
        src = file.split('_')[1]
        bkg = [f for f in os.listdir(self.src_path) if src in f and 'bkg' in f and 'rundata' in f][0]
        return os.path.join(self.src_path,bkg)
    def xray_config(self, file):
        file = os.path.join(self.x_path,file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending='xray', config_file=self.x_config, time_cut=self.time_cut[basename])
        bkg_read_config = file_lib.Read_config(file, ending='xray', config_file=self.x_config, time_cut=self.bkg_time_cut[basename])
        spectrum_config = file_lib.Spectrum_config(corr=self.corr, bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename],self.bkg_form[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]
    
    def src_config(self, file):
        file = os.path.join(self.src_path,file)
        basename = os.path.basename(file)
        bkg = self.__get_src_bkg(basename)
        read_config = file_lib.Read_config(file, ending='normal')
        bkg_read_config = file_lib.Read_config(bkg, ending='normal')
        spectrum_config = file_lib.Spectrum_config(corr=self.corr, bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename],self.bkg_form[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]
    
    def to_x_op(self):
        return Operation(self.x_list, lambda file: self.xray_config(file), self)
    def to_src_op(self):
        return Operation(self.src_list, lambda file: self.src_config(file), self)
    def __center_fit(self, energy:Float1D, center:Float1D, center_err:Float1D):
        popt, pcov = np.polyfit(center, energy,deg=2, full=False, cov=True, w=1./center_err)
        perr = np.sqrt(np.diag(pcov))
        return list(popt), list(perr)
    def __resolution_fit(self, energy, resolution, resolution_err):
        popt, perr = util.resolution_fit(energy, resolution, resolution_err)
        return popt, perr
    def ec_fit(self,src_result,src_energy,x_result,x_energy):
        result = src_result + x_result
        energy = src_energy + x_energy
    
        data = sorted(list(zip(energy, result)), key= lambda x: x[0])
        energy = np.array(list(map(lambda x: x[0],data)))
        center = [np.array([fit[i]['b'] for _,fit in data]) for i in range(4)]
        center_err = [np.array([fit[i]['b_err'] for _,fit in data]) for i in range(4)]
        resolution = [np.array([fit[i]['resolution'] for _,fit in data]) for i in range(4)]
        resolution_err = [np.array([fit[i]['resolution_err'] for _,fit in data]) for i in range(4)]
        q_low = energy<self.energy_split
        q_high = energy>=self.energy_split
        result = [{},{},{},{}]
        for i in range(4):
            ec_low,ec_low_err = self.__center_fit(energy[q_low], center[i][q_low], center_err[i][q_low])
            ec_high,ec_high_err = self.__center_fit(energy[q_high], center[i][q_high], center_err[i][q_high])
            resolution_low, resolution_low_err = self.__resolution_fit(energy[q_low], resolution[i][q_low], resolution_err[i][q_low])
            resolution_high, resolution_high_err = self.__resolution_fit(energy[q_high], resolution[i][q_high], resolution_err[i][q_high])
            
            result[i] = {
                'channel':i,
                'EC_low':ec_low,
                'EC_low_err': ec_low_err,
                'EC_high':ec_high,
                'EC_high_err': ec_high_err,
                'resolution_low':resolution_low,
                'resolution_low_err':resolution_low_err,
                'resolution_high':resolution_high,
                'resolution_high_err':resolution_high_err
            }
            fit_name = f'ec_coef_sci_ch{i}.json'
            util.json_save(result[i], f'{self.result_path}/{util.headtime(fit_name)}')
            save_data = np.array([energy, center[i]], dtype=np.float64)
            data_name = f'ec_data_ch{i}.npy'
            np.save(f'{self.result_path}/{util.headtime(data_name)}', arr=save_data)
        util.ec_plot(energy,center, result, src_energy, x_energy, src_result, x_result, self.result_path)
        return result

        

class TB_operation_03B(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        self.files = [file for file in os.listdir(self.path) if "rundata" in file and 'baseline' not in file and "CI" not in file and '50C' not in file]
        self.adc_max = 16384.0
        self.source = "Am241"
        
        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
    def file_config(self, file):
        file = os.path.join(self.path,file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending='03b')
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]

class EC_operation_03B(EC_operation_05B):
    def __init__(self, tb_result_path:str, fit_range:str, energy:str, bkg_form:str, x_path, src_path, save_path, save_fig_path, result_path) -> None:
        # read config
        self.x_config = ""

        # spectrum config
        ref_temp = 25
        ref_bias = 28.5 
        tb_result:List[Dict[str, Any]] = util.json_load(tb_result_path)
        ref_func = [lambda t,b: util.tempbias2DFunction(t, b, c['G0'],c['k'],c['V0'],c['b'],c['c']) for c in tb_result]
        self.corr = [lambda t,b: f(ref_temp, ref_bias)/f(t,b) for f in ref_func]
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
        self.x_ch = [f for f in os.listdir(self.x_path) if '_rundata' in f and 'CI' not in f]
        
        self.x_list = list(set([f.split('_')[1] for f in self.x_ch]))
        self.x_list.sort()
        # self.src_list = [f for f in os.listdir(self.src_path) if 'rundata' in f and 'bkg' not in f and 'CI' not in f]
        self.src_list = ['src_Na22_20m_10cm_rundata2021-05-05-15-29-56.dat','src_Am241_5m_10cm_rundata2021-05-05-15-15-48.dat','src_Cs137_12m_10cm_rundata2021-05-05-12-12-27.dat']
        self.src_bkg = ['src_bkg_5m_10cm_rundata2021-05-05-15-54-31.dat',   'src_bkg_5m_10cm_rundata2021-05-05-14-55-58.dat',  'src_bkg_5m_10cm_rundata2021-05-05-12-33-57.dat']
        self.energy = util.json_load(energy)
        self.energy_split = 50.2 # keV, absorption edges of Gd
    def __get_x_files(self, energy_name:str):
        return [[os.path.join(self.x_path, f) for f in self.x_ch if f"{energy_name}_ch{i}" in f][0] for i in range(4)]
    def xray_config(self, energy_name:str):
        read_config = [file_lib.Read_config(ch_file, ending='03b', config_file=self.x_config) for ch_file in self.__get_x_files(energy_name)]
        bkg_read_config = read_config[1:4]
        bkg_read_config.append(read_config[0])
        spectrum_config = file_lib.Spectrum_config(corr=self.corr, bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[energy_name],self.bkg_form[energy_name])
        return [read_config, bkg_read_config, spectrum_config, fit_config]
    def __get_src_bkg(self, name:str):
        return self.src_bkg[self.src_list.index(name)]
    def src_config(self, file):
        bkg = os.path.join(self.src_path, self.__get_src_bkg(file))
        file = os.path.join(self.src_path,file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending='03b-src')
        bkg_read_config = file_lib.Read_config(bkg, ending='03b-src')
        spectrum_config = file_lib.Spectrum_config(corr=self.corr, bin_width=self.bin_width)
        fit_config = file_lib.Fit_config(self.fit_range[basename],self.bkg_form[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]

class TB_operation_07B(TB_operation_05B):
    def __init__(self, path, fit_range, save_path, save_fig_path, result_path) -> None:
        self.path = path
        # 
        self.files = [file for file in os.listdir(self.path) if os.path.splitext(file)[1] == '.txt']
        self.adc_max = 65535.0
        self.source = "Am241"
        
        self.fit_range = util.json_load(fit_range)
        self.bin_width = 6
        self.save_path = save_path
        self.save_fig_path = save_fig_path
        self.result_path = result_path
    def file_config(self, file):
        file = os.path.join(self.path,file)
        basename = os.path.basename(file)
        read_config = file_lib.Read_config(file, ending='07')
        bkg_read_config = file_lib.Read_config()
        spectrum_config = file_lib.Spectrum_config(bin_width=self.bin_width, adc_max=self.adc_max)
        fit_config = file_lib.Fit_config(self.fit_range[basename])
        return [read_config, bkg_read_config, spectrum_config, fit_config]
    
def __get_fp05B(config) -> file_lib.File_operation_05b:
    return file_lib.File_operation_05b(config[0].path, *config)
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
def __get_fp03B(config) -> file_lib.File_operation_05b:
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fps = [file_lib.File_operation_05b(read_config[i].path, read_config[i], bkg_read_config[i], spectrum_config, fit_config) for i in range(4)]
    sci = __dict_4ch_reconstruct([fps[i].sci for i in range(4)])
    tel = __dict_4ch_reconstruct([fps[i].tel for i in range(4)])
    bkg_sci = __dict_4ch_reconstruct([fps[i].bkg_sci for i in range(4)])
    bkg_tel = __dict_4ch_reconstruct([fps[i].bkg_tel for i in range(4)])
    fp = fps[0]
    fp.sci, fp.tel = sci, tel
    fp.bkg_sci, fp.bkg_tel = bkg_sci, bkg_tel
    return fp
def process(op:Operation, file:str, fp_method = None,**kw_args) -> None:
    '''process file in sinlge experiment, like specific temp bias or photon energy
    
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
    '''
    fp_method = __get_fp05B if fp_method is None else __get_fp03B
    config = op.file_config(file)
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fp = fp_method(config)
    fp.get_spectrum()
    if kw_args.get("x_lim", None) is not None:
        util.raw_plot(fp.spectrum, fp.x, title=file, x_lim=kw_args["x_lim"], save_path=kw_args.get("save_path", None))
        return 
    fp.peak_fit()
    file = os.path.splitext(file)[0]
    util.fit_plot(fp.spectrum, fp.x, fp.fit_result, title=file, bkgForm=fit_config.bkg_form, fit_range=fit_config.fit_range, save_path=f"{op.op.save_fig_path}/{file}.png")
    fp.save(os.path.join(op.op.save_path,f'{file}.pickle'))