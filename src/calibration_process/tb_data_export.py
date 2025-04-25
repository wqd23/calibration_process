from . import util_lib as util
import numpy as np
from collections import namedtuple, defaultdict
from pathlib import Path


FIT = namedtuple("FIT", ["center", "center_err", "res", "res_err"])

ver = ["03B", "04", "05B", "07"]

path = [Path(f"data/{v}/single_process/TB_fit_result") for v in ver]


def fit_extract(fit_result):
    if fit_result is None:
        return None
    else:
        return FIT(
            fit_result["b"],
            fit_result["b_err"],
            fit_result["resolution"],
            fit_result["resolution_err"],
        )


def extract(path):
    data = util.pickle_load(path)
    fit_result = list(map(fit_extract, data["fit_result"]))
    temp = list(map(np.mean, data["tel"]["tempSipm"]))
    temp_err = list(map(np.std, data["tel"]["tempSipm"]))
    bias = list(map(np.mean, data["tel"]["bias"]))
    bias_err = list(map(np.std, data["tel"]["bias"]))
    file = data["file"]
    result = [
        {
            "fit": fit_result[i],
            "temp": temp[i],
            "temp_err": temp_err[i],
            "bias": bias[i],
            "bias_err": bias_err[i],
        }
        for i in range(4)
    ]
    return file, result


tb_data = defaultdict(lambda: [[], [], [], []])

for v, p in zip(ver, path):
    for file in p.iterdir():
        file, result = extract(file)
        for i in range(4):
            tb_data[v][i].append(
                {
                    "file": file,
                    "fit": result[i]["fit"],
                    "temp": result[i]["temp"],
                    "temp_err": result[i]["temp_err"],
                    "bias": result[i]["bias"],
                    "bias_err": result[i]["bias_err"],
                }
            )
print(tb_data["03B"])
util.pickle_save(tb_data, "./tb_data.pkl")
