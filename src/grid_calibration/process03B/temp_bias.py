import os
from .. import file_lib
from .. import util_lib as util
from ..operation import TB_operation_03B
import argparse
from .__init__ import tb_op
def process(file):
    config = tb_op.file_config(file)
    read_config, bkg_read_config, spectrum_config, fit_config = config
    fp = file_lib.File_operation_05b(config[0].path, *config)
    fp.get_spectrum()
    file = os.path.splitext(os.path.basename(file))[0]
    try:
        fp.peak_fit()
        util.fit_plot(fp.spectrum, fp.x, fp.fit_result, title=file, bkgForm=fit_config.bkg_form, fit_range=fit_config.fit_range, save_path=f"{tb_op.save_fig_path}/{file}.png")
    except Exception as e:
        print(e)
        util.raw_plot(fp.spectrum, fp.x, title=f'{os.path.basename(file)}')
    fp.save(os.path.join(tb_op.save_path,f'{file}.pickle'))
def list_files(args):
    for i, file in enumerate(tb_op.files):
        print(f'{i}: {file}')
def run_single(args):
    n = args.n
    if n == 'all':
        for file in tb_op.files:
            process(file)
        return
    else:
        n = int(n)
        file = tb_op.files[n]
        process(file)
# argparse
parser = argparse.ArgumentParser(description='Process temp bias data')
# add sub command, list, run list_files
subparsers = parser.add_subparsers(help='sub-command help')
list_parser = subparsers.add_parser('list', help='list help')
list_parser.set_defaults(func=list_files)
run_parser = subparsers.add_parser('run', help='run help')
run_parser.add_argument('n', help='一个int参数n或"all"')
run_parser.set_defaults(func=run_single)
if __name__ == "__main__":
    # parser args
    args = parser.parse_args()
    args.func(args)