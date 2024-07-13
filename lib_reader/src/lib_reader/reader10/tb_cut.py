from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def get_tb_cut(utc:np.ndarray, bias:np.ndarray, bias_set:float, tol=0.2):
    # sci collect time
    time = utc[-1] - utc[0]
    # bias set finished, longest should used, but this may be enough
    bias_set_finished = np.logical_and(bias > bias_set - tol, bias < bias_set + tol)
    bias_set_finished = np.where(bias_set_finished)[0]
    t = ((bias_set_finished[-1]-bias_set_finished[0]) - time) // 2
    left, right = bias_set_finished[0] + t, bias_set_finished[0] + t + time
    q = slice(left, right)
    return q

def extract_bias(file:str):
    if 'C_' in file:
        # 温度偏压实验数据
        t = file.split('_')[-1].split('.')[0]
        return int(t)/10.
    else:
        return 28.5

def tel_cut(sci, tel, file:str):
    bias_set = extract_bias(file)
    for ch in range(4):
        utc, bias = sci['utc'][ch], tel['bias'][ch]
        q = get_tb_cut(utc, bias, bias_set)
        tel['bias'][ch] = bias[q]
        tel['tempSipm'][ch] = tel['tempSipm'][ch][q]
    return tel


def plot_tb_cut(sci, tel, file:Path):
    slices = []
    bias_set = extract_bias(str(file))
    fig, axs = plt.subplots(4,1, figsize=(10,10), sharex='all')
    axs[0].set_title(f'{file.stem}: {bias_set}V')
    for ch in range(4):
        plt.sca(axs[ch])
        utc, bias = sci['utc'][ch], tel['bias'][ch]
        q = get_tb_cut(utc, bias, bias_set)
        idx = np.arange(len(bias))
        plt.plot(idx, bias, label='bias/V')
        plt.plot(idx[q], bias[q], 'o', color='red', label='cut')
        plt.axvline(idx[q][0], color='red', linestyle='--')
        plt.axvline(idx[q][-1], color='red', linestyle='--')
        plt.legend()
        slices.append(q)
    fig.tight_layout()
    fig.savefig(f'cut_fig/{file.stem}.png')
    plt.close(fig)
    return slices