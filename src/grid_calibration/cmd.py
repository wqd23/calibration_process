# -*- coding:utf-8 -*-
"""
command line interface for calibration function
----------
"""
from .operation import Operation, process


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
