import numpy as np
from addict import Dict
from pathlib import Path

from ..packet_parser import parse_grid_data_new
from .tb_cut import tel_cut
from ..l1_cache import with_l1_cache, get_l2_processed

_XML = str(Path(__file__).with_name("grid_packet.xml"))


def temp_rebuild(temp:np.ndarray):
    # rebuild raw temp data to unit ℃
    t = temp/16/16
    t = t - (t>128)*255
    assert not ((t>60) | (t<-30)).any(), "temp error. temp > 60 or temp < -30"
    return t


def _readSci_impl(path):
    return Dict(parse_grid_data_new(path, xml_file=_XML, data_tag='grid1x_wf_packet', endian='MSB')[0])


def _readHK_impl(path):
    return Dict(parse_grid_data_new(path, xml_file=_XML, data_tag='hk_grid1x_packet', endian='MSB')[0])


@with_l1_cache(ver="10B", reader="10b", kind="sci")
def readSci(path):
    return _readSci_impl(path)


@with_l1_cache(ver="10B", reader="10b", kind="hk")
def readHK(path):
    return _readHK_impl(path)


def _single_read10_impl(path, hk_name, overwrite):
    wf_data_l = readSci(path, overwrite_cache=overwrite)
    hk_data = readHK(hk_name, overwrite_cache=overwrite)
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

    sciExtracted.pop('waveform_data', None)
    return sciExtracted, telExtracted


def single_read10(path: str, **kwargs):
    overwrite = kwargs.get('overwrite_cache', False)

    def process():
        return _single_read10_impl(path, path.replace('observe', 'hk'), overwrite)

    return get_l2_processed("10B", "10b", path, {}, process, overwrite=overwrite)
