import os
import numpy as np


def get_tb_cut(utc_sci:np.ndarray, utc_tel:np.ndarray):
    q_sci = (utc_tel > utc_sci[0]) & (utc_tel < utc_sci[-1])
    return q_sci


def extract_bias(file:str):
    # 去除idx
    cuts = file.split('_')[1:]
    bias = ["26.5", "27","27.5", "28", "28.3","28.5", "28.7", "29"]
    for b in bias:
        if b in cuts:
            return float(b)
    else:
        raise ValueError(f"bias of {file} not found")

def cut_with_bias(utc:np.ndarray, bias:list[np.ndarray], bias_set:float, tol=0.2):
    # sci collect time
    time = utc[-1] - utc[0]
    bias_set_finished = []
    for b in bias:
        # bias set finished, longest should used, but this may be enough
        finished = np.logical_and(b > bias_set - tol, b < bias_set + tol)
        bias_set_finished.append(finished)
    bias_set_finished = np.logical_and.reduce(bias_set_finished, axis=0)
    bias_set_finished = np.where(bias_set_finished)[0]
    t = ((bias_set_finished[-1]-bias_set_finished[0]) - time) // 2
    left, right = bias_set_finished[0] + t, bias_set_finished[0] + t + time
    q = slice(left, right)
    return q
def tel_cut(sci, tel, file:str):
    utc_sci, utc_tel = sci['utc'], tel['utc_time']
    q = get_tb_cut(utc_sci, utc_tel)
    if not np.any(q):
        bias_set = extract_bias(file)
        q = cut_with_bias(utc_tel, tel['bias'], bias_set)
        print(f"Warning: no utc cut found, using bias cut for {file}")
        raise ValueError(f"no utc cut found for {file}, using bias cut")
    for k in tel.keys():
        if isinstance(tel[k], list) and len(tel[k]) == 4:
            tel[k] = [tel[k][i][q] for i in range(4)]
        elif isinstance(tel[k], np.ndarray):
            tel[k] = tel[k][q]
    return tel