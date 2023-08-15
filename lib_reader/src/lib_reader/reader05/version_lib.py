# -*- coding:utf-8 -*-
"""
GRID version concerining Functions Library
----------

Basic functions for different version of GRID data processing\n
"""
from .gridProcessFunctions import gridProcessFunctions03 as process
from .gridParameters import gridParameters03 as specificParam
from .my_type import *

def refactor_to4chan(v):
    return [np.array(v[i]) for i in range(4)]
def data_refactor(sci:dict):
    return {k:v if len(v) != 4 else refactor_to4chan(v) for k,v in sci.items()}

def single_read05b_normal(path:str)-> Tuple[Dict[str, Float_array_4channel], Dict[str, Float_array_4channel]]:
    """read all file in single temp bias experiment

        Parameters
        ----------
        path : str
            full path of the rundata file

        Returns
        -------
        Tuple[Dict[str, Any],Dict[str, Any]]
            sciExtracted = {
                'amp': array of arrays, extracted amplitude data,
                'timestampEvt': array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
                'eventID': array of arrays, extracted event ID,
                'sciNum': array of arrays, corresponding science run number of each recorded event, for binary files only,
            }
            telExtracted = {
                'tempSipm': array of arrays, extracted SiPM temperature,
                'iMon': array of arrays, extracted monitored current,
                'bias': array of arrays, extracted bias values,
                'iSys': array of arrays, extracted monitored system current,
                'timestamp': array-like or array of arrays, extracted timestamp values,
                'utc': array-like or array of arrays, extracted UTC values, in the same form as 'timestamp',
                'telNum': array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp',
            }
    """
    sciExtracted, telExtracted = process.dataReadout(path, configFile='', singleFile=False, featureMode=False, newDatapack=True, timeCut=None, doDeduplicate=False,
                        useFitBaseline=True, maxUdpReadout=-1, reqNum=[], ignoreCrc=False, plot=False, energyCut=None, baseLen=specificParam.baseLen,
                        noUdp=True, ending="normal")
    sciExtracted, telExtracted = data_refactor(sciExtracted), data_refactor(telExtracted)
    
    return sciExtracted, telExtracted
def single_read05b_xray(path:str, config:str, time_cut=None)-> Tuple[Dict[str, Float_array_4channel], Dict[str, Any]]:
    """read rundata, HK,timeline, config of 05b x_ray style data file
        not pure

    Parameters
    ----------
    path : str
        full path of rundata
    config : str
        config file
    Returns
    -------
    Tuple[Dict[str, Any],Dict[str, Any]]
        sciExtracted = {
            'amp': array of arrays, extracted amplitude data,
            'timestampEvt': array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
            'eventID': array of arrays, extracted event ID,
            'sciNum': array of arrays, corresponding science run number of each recorded event, for binary files only,
        }
        telExtracted = {
            'tempSipm': array of arrays, extracted SiPM temperature,
            'iMon': array of arrays, extracted monitored current,
            'bias': array of arrays, extracted bias values,
            'iSys': array of arrays, extracted monitored system current,
            'timestamp': array-like or array of arrays, extracted timestamp values,
            'utc': array-like or array of arrays, extracted UTC values, in the same form as 'timestamp',
            'telNum': array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp',
        }
    """    
    sciExtracted, telExtracted = process.dataReadout(path, configFile=config, singleFile=False, featureMode=False, newDatapack=True, timeCut=None, doDeduplicate=False,
                        useFitBaseline=True, maxUdpReadout=-1, reqNum=[], ignoreCrc=False, plot=False, energyCut=None, baseLen=specificParam.baseLen,
                        noUdp=True, ending="x_ray")
    sciExtracted, telExtracted = data_refactor(sciExtracted), data_refactor(telExtracted)
    return sciExtracted, telExtracted

def single_read03b(path:str, config:str, time_cut=None):
    sciExtracted, telExtracted = process.dataReadout(path, configFile=config, singleFile=False, featureMode=True, newDatapack=True, timeCut=None, doDeduplicate=False,
                        useFitBaseline=True, maxUdpReadout=specificParam.maxUdpReadout, reqNum=[], ignoreCrc=False, plot=False, energyCut=None, baseLen=specificParam.baseLen,
                        noUdp=False, ending="03b")
    sciExtracted, telExtracted = data_refactor(sciExtracted), data_refactor(telExtracted)
    return sciExtracted, telExtracted
def src_read03b(path:str, config:str, time_cut=None):
    sciExtracted, telExtracted = process.dataReadout(path, configFile=config, singleFile=False, featureMode=False, newDatapack=True, timeCut=None, doDeduplicate=False,
                        useFitBaseline=True, maxUdpReadout=specificParam.maxUdpReadout, reqNum=[], ignoreCrc=False, plot=False, energyCut=None, baseLen=specificParam.baseLen,
                        noUdp=False, ending="03b")
    sciExtracted, telExtracted = data_refactor(sciExtracted), data_refactor(telExtracted)
    return sciExtracted, telExtracted