from .parse_grid_data import parse_grid_data_new
from addict import Dict
from pathlib import Path
from ..reader10.read import temp_rebuild
from ..reader10.tb_cut import tel_cut

def readSci(path, mode='wf'):
    if mode == 'wf':
        return Dict(parse_grid_data_new(path,data_tag='grid1x_wf_packet',endian='MSB')[0])
    elif mode == 'hk':
        return Dict(parse_grid_data_new(path,data_tag='grid1x_ft_packet',endian='MSB')[0],multi_evt=41, multi_step=12)
    else:
        raise ValueError("Invalid mode. Use 'wf' or 'hk'.")

def readHK(path):
    return Dict(parse_grid_data_new(path,data_tag='hk_grid1x_packet',endian='MSB')[0])

def getHK(sciFile):
    sciFile = Path(sciFile)
    hkFile1 = sciFile.with_name(sciFile.name.replace('observe', 'hk'))
    hkFile2 = sciFile.with_name(sciFile.name.replace('_observe.dat', '.hk'))
    hkFile3 = sciFile.with_name(sciFile.name.replace('_observe.dat', 'hk'))
    if not any([hkFile1.exists(), hkFile2.exists(), hkFile3.exists()]):
        raise FileNotFoundError(f"HK file not found for {sciFile}.")
    return hkFile1 if hkFile1.exists() else hkFile2 if hkFile2.exists() else hkFile3

def single_read11(path:str, mode='wf'):
    observe_name = path
    hk_name = getHK(observe_name)
    wf_data_l = readSci(observe_name, mode=mode)
    hk_data = readHK(hk_name)
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
    telExtracted.tempSipm = [telExtracted[f'sipm_temp{i}']/100-273.15 for i in range(4)]
    # current, unit uA
    telExtracted.iMon = [telExtracted[f'sipm_current{i}'] for i in range(4)]
    # bias monitor, unit V
    telExtracted.vMon = [telExtracted[f'sipm_voltage{i}']/1000 for i in range(4)]
    telExtracted.bias = [telExtracted.vMon[i] - 499*telExtracted.iMon[i]*1e-6  for i in range(4)]
    
    sciExtracted['timestampEvt'] = sciExtracted.timestamp

    # time cut
    telExtracted = tel_cut(sciExtracted, telExtracted, path)

    return sciExtracted, telExtracted