from .parse_grid_data import parse_grid_data_new
from addict import Dict
from pathlib import Path
from ..reader11.tb_cut import tel_cut
from cachier import cachier

@cachier(cache_dir=Path('.cache') / '11B', separate_files=True)
def readSci(path, mode='wf'):
    if mode == 'wf':
        return Dict(parse_grid_data_new(path,data_tag='grid1x_wf_packet',endian='MSB')[0])
    elif mode == 'ft':
        return Dict(parse_grid_data_new(path,data_tag='grid1x_ft_packet',endian='MSB')[0],multi_evt=41, multi_step=12)
    else:
        raise ValueError("Invalid mode. Use 'wf' or 'hk'.")

def readHK(path):
    return Dict(parse_grid_data_new(path,data_tag='hk_grid1x_packet',endian='MSB')[0])

def getHK(sciFile):
    sciFile = Path(sciFile)
    idx = sciFile.stem.split('_')[0]
    # temp-bias, src data
    Files = [f for f in sciFile.parent.glob(f'*')]
    hkFile = [f for f in Files if f.stem.startswith(idx) and 'hk' in f.name]
    if len(hkFile) == 0:
        raise FileNotFoundError(f"HK file {hkFile} does not exist.")
    hkFile = hkFile[0]
    return hkFile

def single_read11(path:str, mode='wf', **kwargs):
    observe_name = path
    hk_name = getHK(observe_name)
    wf_data_l = readSci(observe_name, mode=mode, overwrite_cache=kwargs.get('overwrite_cache', False))
    hk_data = readHK(hk_name)
    sciExtracted, telExtracted = wf_data_l, hk_data
    
    # amp
    if len(sciExtracted.data_max) == len(sciExtracted.data_base):
        amp = sciExtracted.data_max - sciExtracted.data_base/4.
    else:
        assert False, f'{path} data_max.len != data_base.len'
    sciExtracted.amp = amp

    telExtracted.tempSipm = [telExtracted[f'sipm_temp{i}']/100-273.15 for i in range(4)]
    # current, unit uA
    telExtracted.iMon = [telExtracted[f'sipm_current{i}'] for i in range(4)]
    # bias monitor, unit V
    telExtracted.vMon = [telExtracted[f'sipm_voltage{i}']/1000 for i in range(4)]
    telExtracted.bias = [telExtracted.vMon[i] - 499*telExtracted.iMon[i]*1e-6  for i in range(4)]
    # tempSipm
    
    sciExtracted['timestampEvt'] = sciExtracted.timestamp

    # time cut
    telExtracted = tel_cut(sciExtracted, telExtracted, path)

    n = sciExtracted.data_max.shape[0]
    for k,v in sciExtracted.items():
        if v.shape[0] == n and k != 'channel_n':
            sciExtracted[k] = [v[sciExtracted.channel_n == i] for i in range(4)]
    del n

    return sciExtracted, telExtracted