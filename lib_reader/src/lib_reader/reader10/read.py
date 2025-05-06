import numpy as np
from addict import Dict

from ..util import data_refactor
from .parse_grid_data import parse_grid_data_new
from .tb_cut import tel_cut


def temp_rebuild(temp:np.ndarray):
    # rebuild raw temp data to unit ℃
    t = temp/16/16
    t = t - (t>128)*255
    assert not ((t>60) | (t<-30)).any(), "temp error. temp > 60 or temp < -30"
    return t

def single_read10(path:str):
    observe_name = path
    hk_name = observe_name.replace('observe', 'hk')
    wf_data_l = Dict(parse_grid_data_new(observe_name,data_tag='grid1x_wf_packet',endian='MSB')[0])
    hk_data = Dict(parse_grid_data_new(hk_name,data_tag='hk_grid1x_packet',endian='MSB')[0])
    sciExtracted, telExtracted = wf_data_l, hk_data
    
    # amp
    if len(sciExtracted.data_max) == len(sciExtracted.data_base):
        amp = sciExtracted.data_max - sciExtracted.data_base/4.
    else:
        assert False, f'{path} data_max.len != data_base.len'
    sciExtracted.amp = amp

    n = sciExtracted.data_max.shape[0]
    for k,v in sciExtracted.items():
        if v.shape[0] == n and k != 'channel_n':
            sciExtracted[k] = [v[sciExtracted.channel_n == i] for i in range(4)]
    del n

    # tempSipm
    telExtracted.tempSipm = [temp_rebuild(telExtracted[f'sipm_temp{i}']) for i in range(4)]
    # current, unit uA
    telExtracted.iMon = [2.5*telExtracted[f'sipm_current{i}']/4096/548.88*1000_000 for i in range(4)]
    # bias monitor, unit V
    telExtracted.vMon = [20.57*2.5*telExtracted[f'sipm_voltage{i}']/4096 for i in range(4)]
    telExtracted.bias = [20.57*2.5*telExtracted[f'sipm_voltage{i}']/4096 - 499*2.5*telExtracted[f'sipm_current{i}']/4096/548.88 for i in range(4)]
    
    sciExtracted['timestampEvt'] = sciExtracted.timestamp

    # time cut
    telExtracted = tel_cut(sciExtracted, telExtracted, path)

    return sciExtracted, telExtracted