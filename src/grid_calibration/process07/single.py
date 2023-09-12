from lib_plot import plot
from ..import file_lib
from .. import util_lib as util
import os
path = "data/07/raw_data/北师大正样温度偏压标定数据-20211124/211120150442_COM3-bnu_tb_26p5V_0C_Thres6_90s.txt"
reader_config = file_lib.Read_config(path, ending='07')
bkg_reader_config = file_lib.Read_config()
spectrum_config = file_lib.Spectrum_config()
fit_config = file_lib.Fit_config([[400,1000],[400,1000],[500,1000],[400,800]])
fp = file_lib.File_operation_05b(path,reader_config, bkg_reader_config, spectrum_config, fit_config)
fp.get_spectrum()
plot.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(path)}', save_path='raw.png', x_lim = (0,3000))
fp.peak_fit()
plot.fit_plot(fp.spectrum, fp.x, fp.fit_result, title=f'{os.path.basename(path)}',bkgForm='lin', save_path='fit.png')