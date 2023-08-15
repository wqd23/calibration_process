from ..import file_lib
from .. import util_lib as util
import os
path = "data/03B/raw_data/20210501_tempbias_03B/-10C_265_4m_rundata2021-05-01-14-27-37.dat"
reader_config = file_lib.Read_config(path, ending='03b')
bkg_reader_config = file_lib.Read_config()
spectrum_config = file_lib.Spectrum_config()
fit_config = file_lib.Fit_config([None, None, [197, 393], [221, 468]])
fp = file_lib.File_operation_05b(path,reader_config, bkg_reader_config, spectrum_config, fit_config)
fp.get_spectrum()
util.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(path)}', save_path='raw.png')
fp.peak_fit()
util.fit_plot(fp.spectrum, fp.x, fp.fit_result, title=f'{os.path.basename(path)}',bkgForm='lin', save_path='fit.png')
fp.save('.')