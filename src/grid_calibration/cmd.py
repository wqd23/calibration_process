# -*- coding:utf-8 -*-
"""
command line interface for calibration function
----------
"""
from .operation import Operation, process
from . import util_lib as util
class file_list_op:
    def __init__(self, op: Operation, fp_method=None) -> None:
        self.files = op.files
        self.file_config = op.file_config
        self.op = op
        self.fp_method = fp_method

    def __get_file(self, n):
        try:
            file = self.files[n]
        except IndexError as e:
            print(f"length of file list if {len(self.files)}, but {n} got")
            raise e
        return file

    # list all the file could be processed
    def list(self):
        for i, file in enumerate(self.files):
            print(i, file)

    # plot the raw figure of {n}th file, within(left, right)
    def raw(self, n, left, right, save_path=None):
        process(
            self.op,
            self.__get_file(n),
            fp_method=self.fp_method,
            x_lim=[left, right],
            save_path=save_path,
        )

    def __process_list(self, n_list: list):
        for n in n_list:
            process(self.op, self.__get_file(n), fp_method=self.fp_method)

    # process the {n}th file of list, or "all" file
    def run(self, n, end=None):
        if n == "all":
            self.__process_list(range(len(self.files)))
        else:
            if end is None:
                self.__process_list([n])
            else:
                self.__process_list(range(n, end + 1))

class VersionProcessOp:
    '''
    commandline warpper for different payload
    '''
    def __init__(self, tb_op, ec_op, fp_method) -> None:
        self.tb = file_list_op(tb_op.to_op())
        self.ec = {
            "x": file_list_op(ec_op.to_x_op(), fp_method=fp_method), 
            "src":file_list_op(ec_op.to_src_op())}
        self.__tb_op = tb_op
        self.__ec_op = ec_op
    def tbfit(self):
        data_all = self.__tb_op.load_data()
        self.__tb_op.temp_bias_fit(data_all)
    def ecfit(self):
        x_res, src_res = util.get_fit_dict(self.__ec_op.save_path, self.__ec_op.energy)

        src_result = list(src_res.values())
        src_energy = list(src_res.keys())
        x_result = list(x_res.values())
        x_energy = list(x_res.keys())
        self.__ec_op.ec_fit(src_result, src_energy, x_result, x_energy)