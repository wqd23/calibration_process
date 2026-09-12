from ..packet_parser import parse_grid_data_new
from addict import Dict
from pathlib import Path
from ..reader11.tb_cut import tel_cut
from ..l1_cache import with_l1_cache, get_l2_processed, apply_selection

_XML = str(Path(__file__).with_name("grid_packet.xml"))

def _readSci_impl(path, mode='wf'):
    if mode == 'wf':
        return Dict(parse_grid_data_new(path,xml_file=_XML,data_tag='grid1x_wf_packet',endian='MSB')[0])
    elif mode == 'ft':
        return Dict(parse_grid_data_new(path,xml_file=_XML,data_tag='grid1x_ft_packet',endian='MSB')[0],multi_evt=41, multi_step=12)
    else:
        raise ValueError("Invalid mode. Use 'wf' or 'hk'.")

def _readHK_impl(path):
    return Dict(parse_grid_data_new(path,xml_file=_XML,data_tag='hk_grid1x_packet',endian='MSB')[0])

@with_l1_cache(ver="11B", reader="11b", kind="sci")
def readSci(path, mode='wf'):
    return _readSci_impl(path, mode=mode)

@with_l1_cache(ver="11B", reader="11b", kind="hk")
def readHK(path):
    return _readHK_impl(path)

def getHK(sciFile):
    sciFile = Path(sciFile)
    idx = sciFile.stem.split('_')[0]
    # temp-bias, src data
    Files = [f for f in sciFile.parent.glob('*')]
    hkFile = [f for f in Files if f.stem.startswith(idx) and 'hk' in f.name]
    if len(hkFile) == 0:
        raise FileNotFoundError(f"HK file {hkFile} does not exist.")
    hkFile = hkFile[0]
    return hkFile

def _single_read11_impl(path, mode, hk_name, overwrite, select=None):
    wf_data_l = readSci(path, mode=mode, overwrite_cache=overwrite)
    if select is not None:
        wf_data_l = apply_selection("11B", "11b", path, {"mode": mode}, select[0], select[1], wf_data_l)
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

    sciExtracted.pop('waveform_data', None)
    return sciExtracted, telExtracted


def single_read11(path: str, config=None, mode='wf', **kwargs):
    overwrite = kwargs.get('overwrite_cache', False)
    select = kwargs.get('select')
    l2_params = {"mode": mode}
    if select is not None:
        l2_params["select"] = select[0]

    def process():
        return _single_read11_impl(path, mode, getHK(path), overwrite, select=select)

    return get_l2_processed("11B", "11b", path, l2_params, process,
                            overwrite=overwrite)