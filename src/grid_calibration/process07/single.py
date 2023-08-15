from ..import file_lib
from .. import util_lib as util
import os
path = "data/07/raw_data/北师大豁免源标定数据-20211127/211127154145_COM3-bnu_src_Cs137_45min_10cm.txt"
reader_config = file_lib.Read_config(path, ending='07')
bkg_reader_config = file_lib.Read_config()
spectrum_config = file_lib.Spectrum_config()
fit_config = file_lib.Fit_config([0,10]*4)
fp = file_lib.File_operation_05b(path,reader_config, bkg_reader_config, spectrum_config, fit_config)
fp.get_spectrum()
util.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(path)}', save_path='test.png')
