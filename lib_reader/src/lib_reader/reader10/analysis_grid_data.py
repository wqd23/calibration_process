import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from datetime import datetime

## usage
# datafile = r'/home/liping/grid_data/experiments/grid06b/20230306/200138482_JL1MF02A04_tiange/file-5_IV.bin'
# data_iv,_ = [Dict(v) for v in parse_grid_data_new(datafile,endian='LSB',data_tag='iv_packet')]
# data_vbr,_ = [Dict(v) for v in parse_grid_data_new(datafile,endian='LSB',data_tag='vbr_packet')]
# plot_iv_vbr_curve(cal_sipm_iv_vbr_data(data_iv.iv,scan_type='iv'),cal_sipm_iv_vbr_data(data_vbr.vbr,scan_type='vbr'))

def cal_sipm_iv_vbr_data(data_raw, scan_type='iv', points=25, channel=4):
    '''
    unit: V, $\mu A$
    '''
    Vval = data_raw.reshape(points,channel,2)[...,0]
    Ival = data_raw.reshape(points,channel,2)[...,1]

    Rval = 0
    if scan_type=='iv':
        Rval = 49.9*2e5/(2e5+49.9)
    elif scan_type=='vbr':
        Rval = 2e5
    else:
        ValueError('undefine scan type')

    current = 2.5*Ival/4096./(499+Rval)
    voltage = 20.57*2.5*Vval/4096 - current*499

    return (voltage,current)

def plot_iv_vbr_curve(data_iv, data_vbr):
    '''
    data shape: (point_number, channel_number, 2(V+I))
    '''
    shape_iv = data_iv[0].shape
    shape_vbr = data_vbr[0].shape

    fig = plt.figure('iv_vbr',clear=True,figsize=(18,6))
    ax = fig.subplots(1,2)
    for ich in range(shape_iv[1]):
        ax[0].plot(data_iv[0][:22,ich], data_iv[1][:22,ich]*1e6, '.', label=f'ch-{ich}')
    ax[0].set_title('IV Scan')
    for ich in range(4):
        ax[1].plot(data_vbr[0][:22,ich], data_vbr[1][:22,ich]*1e6, '.', label=f'ch-{ich}')
    ax[1].set_title('Vbr Scan')

    for a in ax:
        a.set_xlabel('bias voltage (V)')
        a.set_ylabel('leakage current ($\mu A$)')
        a.grid(which='major',ls='--')
        a.legend()
    fig.show()

def cg_time_conv(dt_year,dt_month,dt_day,dt_hour,dt_min=0,dt_sec=0):
    """
    This function converts the cg date and time to Calendar time
    """
    return datetime.fromtimestamp(datetime(dt_year,dt_month,dt_day,dt_hour,dt_min,dt_sec).timestamp()+datetime(2000,1,1,20).timestamp())