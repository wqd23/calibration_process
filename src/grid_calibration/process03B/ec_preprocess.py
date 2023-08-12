from .. import util_lib as util
import os
from ..operation import EC_operation_03B
from .. import file_lib
import argparse
from .__init__ import ec_op

def dict_4ch_reconstruct(dict_4ch):
    # dict_4ch is list of 4 dict all value to be list with 4 elements, return a dict with 4 all value to be 4 ch 
    dict_4ch_re = {}
    for key in dict_4ch[0].keys():
        dict_4ch_re[key] = []
        for i in range(4):
            dict_4ch_re[key].append(dict_4ch[i][key][i])
    return dict_4ch_re 
def get_fp(config) -> file_lib.File_operation_05b:
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fps = [file_lib.File_operation_05b(read_config[i].path, read_config[i], bkg_read_config[i], spectrum_config, fit_config) for i in range(4)]
    sci = dict_4ch_reconstruct([fps[i].sci for i in range(4)])
    tel = dict_4ch_reconstruct([fps[i].tel for i in range(4)])
    bkg_sci = dict_4ch_reconstruct([fps[i].bkg_sci for i in range(4)])
    bkg_tel = dict_4ch_reconstruct([fps[i].bkg_tel for i in range(4)])
    fp = fps[0]
    fp.sci, fp.tel = sci, tel
    fp.bkg_sci, fp.bkg_tel = bkg_sci, bkg_tel
    return fp
def process_xfile(x_energy):
    config = ec_op.xray_config(x_energy)
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fp = get_fp(config)
    fp.get_spectrum()
    # util.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(x_energy)}')
    try:
        fp.peak_fit()
        # util.fit_plot_single_channel(fp.spectrum[0], fp.x[0], fp.fit_result[0], title=x_file, bkgForm=fit_config.bkg_form, fit_range=fit_config.fit_range[0])
        util.fit_plot(fp.spectrum, fp.x, fp.fit_result, title=x_energy, bkgForm=fit_config.bkg_form, fit_range=fit_config.fit_range, save_path=f'{ec_op.save_fig_path}/{x_energy}.png')
    except Exception as e:
        print(e)
        util.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(x_energy)}')
    data = {
        "file": x_energy,
        "fit_result": fp.fit_result,
        "spectrum": fp.spectrum,
        "x": fp.x,
        "tel": fp.tel
    }
    util.pickle_save(data,f'{ec_op.save_path}/{x_energy}.pickle')
def process_src(src_file):
    config = ec_op.src_config(src_file)
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fp = file_lib.File_operation_05b(config[0].path, *config)
    fp.get_spectrum()
    src_file = os.path.splitext(os.path.basename(src_file))[0]
    # util.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(src_file)}', save_path=f'03B/single_process/raw/{src_file}.png')
    # return
    try:
        fp.peak_fit()
        util.fit_plot(fp.spectrum, fp.x, fp.fit_result, title=src_file, bkgForm=fit_config.bkg_form, fit_range=fit_config.fit_range, save_path=f'{ec_op.save_fig_path}/{src_file}.png')
    except Exception as e:
        print(e)
        util.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(src_file)}')
    data = {
        "file": src_file,
        "fit_result": fp.fit_result,
        "spectrum": fp.spectrum,
        "x": fp.x,
        "tel": fp.tel
    }
    util.pickle_save(data,f'{ec_op.save_path}/{src_file}.pickle')
def file_list(data_type):
    if data_type == 'x':
        file_list = ec_op.x_list
        process_file = process_xfile
    elif data_type == 'src':
        file_list = ec_op.src_list
        process_file = process_src
    else:
        raise ValueError('data type error')
    for i, file in enumerate(file_list):
        print(f"{i} {file}")
def run_xfile(n):
    if n == 'all':
        for i in range(len(ec_op.x_list)):
            process_xfile(ec_op.x_list[i])
    else:
        n = int(n)
        process_xfile(ec_op.x_list[n])
def run_src(n):
    if n == 'all':
        for i in range(len(ec_op.src_list)):
            process_src(ec_op.src_list[i])
    else:
        n = int(n)
        process_src(ec_op.src_list[n])
parser = argparse.ArgumentParser(description='Process ec data')
# add sub command, list, run list_files
subparsers = parser.add_subparsers(help='sub-command help')
xray_parser = subparsers.add_parser('x', help='xray help')
x_subparsers = xray_parser.add_subparsers(help='xray sub-command help')
xray_list_parser = x_subparsers.add_parser('list', help='list xray files')
xray_list_parser.set_defaults(func=lambda _: file_list('x'))
xray_run_parser = x_subparsers.add_parser('run', help='run xray files')
xray_run_parser.add_argument('n', help='x run 子命令的参数')
xray_run_parser.set_defaults(func=lambda x: run_xfile(x.n))

src_parser = subparsers.add_parser('src', help='src help')
src_subparsers = src_parser.add_subparsers(help='src sub-command help')
src_list_parser = src_subparsers.add_parser('list', help='list src files')
src_list_parser.set_defaults(func=lambda _: file_list('src'))
src_run_parser = src_subparsers.add_parser('run', help='run src files')
src_run_parser.add_argument('n', help='src run 子命令的参数')
src_run_parser.set_defaults(func=lambda x: run_src(x.n))
if __name__ == "__main__":
    # parser args
    args = parser.parse_args()
    args.func(args)