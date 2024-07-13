# -*- coding:utf-8 -*-

"""
GRID03 I/O and Basic Process Functions Library
----------
Basic data readout/output and process functions\n
`v1.0.0` GRID03 and GRID03-based detector calibration result analysis, also quick view of GRID02 in-orbit data, by ghz, cjr and ydx\n
"""

#**************************************************************************************************************************************************
#***************************************************************Import packages***************************************************************
#**************************************************************************************************************************************************

from .. import gridParametersCommon as parameters
from .. gridParameters import gridParameters03 as specificParam
from .. import gridBasicFunctions as basic
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import struct
import os
import re
import json, csv
from copy import copy

#****************************************************************************************************************************************************************
#*******************************************************Basic data readout and process functions*******************************************************
#****************************************************************************************************************************************************************

#****************************************************************************************************************************************************************
#***************************************************************EC config readout functions***************************************************************
#****************************************************************************************************************************************************************

def readECconfig(configFile, dataType = ''):

    """
    Function to read EC configurations from specified config file

    Parameters
    ----------
    configFile : str,
        name of the config file (with full path), namely `xrayConfig`, `sourceConfig` and `HPGeConfig` in gridParameters03.py\n
    dataType : str,
        indicates type of data, with 'xray' being x-ray data, 'source' being source data and `hpge` being HPGe data

    Returns
    ----------
    fileList : list of dicts or dict,
        list of files to readout for all 4 channels for x-ray data or dict for source data and HPGe data, with list of each channel being a dictionary \
            in the form: \n
        files = {
                'date' :             list of experiment dates, for x-ray and HPGe data only
                'energy' :          list of corresponding energy, for x-ray and HPGe data only
                'source' :           list of corresponding sources, for source data only
                'data' :              list of data file names,
                'cfg' :                list of corresponding science config files, for x-ray and source data only
                'bkg' :               list of corresponding background data files, for x-ray and source data only
                'bkg_cfg' :         list of corresponding background science config files, for x-ray and source data only
            }

    Raises
    ----------
    Exception
        when given data type is not an available type, or error reading config file (e.g. config file with wrong format)
    """

    typesAvailable = ['xray', 'source', 'hpge']
    if not dataType in typesAvailable:
        raise Exception('readECconfig: data type \"' + dataType + '\" is not an available type')

    config = []
    configDict = []
    for ich in range(4):
        configDict.append({})
    try:
        with open(configFile, 'r', encoding = 'gbk') as fin:
            reader = csv.reader(fin)
            for line in reader:
                config.append(line)

        config = config[1:]

        #X-ray EC readout config
        if dataType == 'xray':
            channelSet = np.where(np.array(config)[:, 0] != '')[0]
            if len(channelSet) > 4 or len(channelSet) == 0:
                raise Exception()
            for ich in range(len(channelSet)):
                channel = int(np.array(config)[:, 1][channelSet[ich]])
                if configDict[channel]:
                    #Duplicated channel
                    raise Exception()
                nextSet = len(config) if ich >= len(channelSet) - 1 else channelSet[ich + 1]
                subSet = np.where(np.array(config)[channelSet[ich]:nextSet, 2] != '')[0]
                date = []
                energy = []
                dataFile = []
                bkgFile = []
                cfgFile = []
                bkgCfgFile = []
                filePath = config[channelSet[ich]][0]
                for iset in range(len(subSet)):
                    fileSetPath = filePath + config[channelSet[ich] + subSet[iset]][2]
                    nextSubSet = nextSet if iset >= len(subSet) - 1 else channelSet[ich] + subSet[iset + 1]
                    currDate = config[channelSet[ich] + subSet[iset]][3]
                    for ifile in range(channelSet[ich] + subSet[iset], nextSubSet):
                        date.append(currDate)
                        energy.append(float(config[ifile][4]))
                        dataFile.append(fileSetPath + config[ifile][5])
                        bkgFile.append(fileSetPath + config[ifile][6])
                        cfgFile.append((filePath + config[ifile][7]) if config[ifile][7] != '' else '')
                        bkgCfgFile.append((filePath + config[ifile][8]) if config[ifile][8] != '' else '')
                configDict[channel] = {
                    'date' :             date, 
                    'energy' :          energy,
                    'data' :             dataFile,
                    'bkg' :               bkgFile,
                    'cfg' :               cfgFile, 
                    'bkg_cfg' :        bkgCfgFile, 
                }

                
        #Source EC readout config
        elif dataType == 'source':
            source = []
            dataFile = []
            bkgFile = []
            cfgFile = []
            bkgCfgFile = []
            for ifile in range(len(config)):
                filePath = config[ifile][0]
                source.append(config[ifile][1])
                dataFile.append(filePath + config[ifile][2])
                bkgFile.append(filePath + config[ifile][3])
                cfgFile.append((filePath + config[ifile][4]) if config[ifile][4] != '' else '')
                bkgCfgFile.append((filePath + config[ifile][5]) if config[ifile][5] != '' else '')
            configDict = {
                'source' :          source,
                'data' :             dataFile,
                'bkg' :               bkgFile,
                'cfg' :               cfgFile, 
                'bkg_cfg' :        bkgCfgFile, 
            }

        #HPGe EC readout config
        else:
            date = []
            energy = []
            dataFile = []
            rootSet = np.where(np.array(config)[:, 0] != '')[0]
            if len(rootSet) == 0:
                raise Exception()
            for iroot in range(len(rootSet)):
                nextSet = len(config) if iroot >= len(rootSet) - 1 else rootSet[iroot + 1]
                subSet = np.where(np.array(config)[rootSet[iroot]:nextSet, 2] != '')[0]
                filePath = config[rootSet[iroot]][0]
                for iset in range(len(subSet)):
                    fileSetPath = filePath + config[rootSet[iroot] + subSet[iset]][1]
                    nextSubSet = nextSet if iset >= len(subSet) - 1 else rootSet[iroot] + subSet[iset + 1]
                    currDate = config[rootSet[iroot] + subSet[iset]][2]
                    for ifile in range(rootSet[iroot] + subSet[iset], nextSubSet):
                        date.append(currDate)
                        energy.append(float(config[ifile][3]))
                        dataFile.append(fileSetPath + config[ifile][4])
            configDict = {
                'date' :             date, 
                'energy' :          energy,
                'data' :             dataFile,
            }
    except:
        raise Exception('readECconfig: error reading config file \"' + configFile + '\"')

    return configDict

def readECProcessConfig(configFile):

    """
    Function to read EC data process configurations from specified config file

    Parameters
    ----------
    configFile : str,
        name of the config file (with full path), namely `xrayProcessConfig` and `sourceProcessConfig` in gridParameters03.py\n
    
    Returns
    ----------
    processOptions : tuple, 
        in the form of `(newDatapack, singleFile, featureMode, timeCut, rateStyle, deduplicate, useFitBaseline, maxUdpReadout, doCorr, odr, maxiter, \
            bound)`\n
    the options include: \n
        newDatapack : boolean, optional
            indicates whether the data is produced with new datapack\n
        singleFile : boolean, optional
            indicates whether the data is stored in seperate files (rundata, HK, timeline) or single file\n
        featureMode : boolean, optional
            indicates whether the input file is obtained with feature mode\n
        timeCut : int
            cut of time data in seconds, specially designed for temp-bias data with pid bias control(6th ver.)\n
        rateStyle : str
            the style of count rate correction, '' for no corraction, 's' for correction with small data packs(512byte), 'l' for calculation \
                with large data packs(4096byte)\n
        deduplicate : boolean
            indicates whether the data will be deduplicated, for BINARY(.dat) data of new programme (6th ver.) only\n
        useFitBaseline : boolean, optional
            indicates whether the amplitude data will be calculated with measured mean baseline, if False the baseline value in the config file will be used\n
        maxUdpReadout : int, optional
            maximum number of UDP packages to readout in a single readout run, for data readout in multiple readout runs in order to process large \
                data packs, if `maxUdpReadout` <= 0 then the data will be readout in a single run (default readout process)\n
        doCorr : boolean
            indicates whether the temperature-bias correction will be done, to avoid warning info in the output\n
        odr : boolean
            indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
        maxiter : int
            maximum number of iterations, 0 for auto range correction off\n
        bound : float
            boundary for auto rangecorrection in `\sigma`, with upper and lower bounds being `\mu - bound * \sigma` and `\mu + bound * \sigma`\n
        noudp : bool
            whether to use udp unpack
    Raises
    ----------
    Exception
        when error reading config file (e.g. config file with wrong format)
    """

    newDatapack = False
    singleFile = False
    featureMode = False
    timeCut = 0.
    rateStyle = ''
    deduplicate = False
    useFitBaseline = False
    maxUdpReadout = -1
    doCorr = False
    odr = False
    maxiter = 1
    bound = 3.
    noudp = False
    try:
        with open(configFile, 'r') as fin:
            config = json.load(fin)
        #Determine whether the individual options are in correct type
        if isinstance(config['newDatapack'], bool):
            newDatapack = config['newDatapack']
        else:
            raise Exception()
        if isinstance(config['singleFile'], bool):
            singleFile = config['singleFile']
        else:
            raise Exception()
        if isinstance(config['featureMode'], bool):
            featureMode = config['featureMode']
        else:
            raise Exception()
        if isinstance(config['timeCut'], float) or isinstance(config['timeCut'], int):
            timeCut = config['timeCut']
        else:
            raise Exception()
        if isinstance(config['rateStyle'], str):
            rateStyle = config['rateStyle']
        else:
            raise Exception()
        if isinstance(config['deduplicate'], bool):
            deduplicate = config['deduplicate']
        else:
            raise Exception()
        if isinstance(config['useFitBaseline'], bool):
            useFitBaseline = config['useFitBaseline']
        else:
            raise Exception()
        if isinstance(config['maxUdpReadout'], int):
            maxUdpReadout = config['maxUdpReadout']
        else:
            raise Exception()
        if isinstance(config['doCorr'], bool):
            doCorr = config['doCorr']
        else:
            raise Exception()
        if isinstance(config['odr'], bool):
            odr = config['odr']
        else:
            raise Exception()
        if isinstance(config['maxiter'], int):
            maxiter = config['maxiter']
        else:
            raise Exception()
        if isinstance(config['bound'], float) or isinstance(config['bound'], int):
            bound = config['bound']
        else:
            raise Exception()
        if isinstance(config['noudp'], bool):
            noudp = config['noudp']
        else:
            raise Exception()
    except:
        raise Exception('readECProcessConfig: error reading config file \"' + configFile + '\"')

    return newDatapack, singleFile, featureMode, timeCut, rateStyle, deduplicate, useFitBaseline, maxUdpReadout, doCorr, odr, maxiter, bound, noudp

#****************************************************************************************************************************************************************
#*******************************************************************Data readout functions*****************************************************************
#****************************************************************************************************************************************************************

def readSciConfig(filename, noUdp = False):

    """
    Function to readout sicence config dictionary

    Parameters
    ----------
    filename : str, 
        name of the config file with full path. Please note that only json file is supported
    noUdp : boolean, optional
        indicates whether the science raw data is arranged in UDP-level packages\n
    
    Returns
    ----------
    configDict : dict, 
        science config read from the file, in the form of a dictionary as: \n
        for old version datapack with UDP packages (noUdp = False): \n
        configDict = {
                'channel0Bias' :                            set bias of channel 0, 
                'channel1Bias' :                            set bias of channel 1, 
                'channel2Bias' :                            set bias of channel 2, 
                'channel3Bias' :                            set bias of channel 3, 
                'channel0TriggerThreshold' :         trigger threshold of channel 0, 
                'channel1TriggerThreshold' :         trigger threshold of channel 1, 
                'channel2TriggerThreshold' :         trigger threshold of channel 2, 
                'channel3TriggerThreshold' :         trigger threshold of channel 3, 
                'channel0BaseLine' :                     set baseline of channel 0, 
                'channel1BaseLine' :                     set baseline of channel 1, 
                'channel2BaseLine' :                     set baseline of channel 2, 
                'channel3BaseLine' :                     set baseline of channel 3, 
                'channel0Vmin' :                           minimum allowed voltage of channel 0, 
                'channel1Vmin' :                           minimum allowed voltage of channel 1, 
                'channel2Vmin' :                           minimum allowed voltage of channel 2, 
                'channel3Vmin' :                           minimum allowed voltage of channel 3, 
                'SamplingDepth' :                          sampling depth of pulse shape, for pulse shape mode only
                'BaselineLenth' :                           baseline measurement length of pulse shape, for pulse shape mode only
            }
        for new version datapack without UDP packages (noUdp = True): \n
        configDict = {
                'bias_value' :                                set bias of all 4 channels in a list, 
                'time_line_cycle' :                        indicates whether the data is acquired in timeline cycle, 
                'detection_cycle' :                        indicates whether the data is acquired in detection cycle, 
                'trigger_threshold' :                      trigger threshold of all 4 channels in a list, 
                'sampling_depth' :                        sampling depth of pulse shape, for pulse shape mode only, 
                'base_line_length' :                      baseline measurement length of pulse shape, for pulse shape mode only, 
                'ref_value' :                                 bias reference values of all 4 channels in a list, 
                'polarity' :                                    signal polarity, 
                'charge_injection_fre' :                frequency for charge injection signals, 
                'charge_injection_width' :            width for charge injection signals, 
                'charge_injection_control' :          indicates whether the charge injection is on or off, 
                'data_ack' :                                  indicates whether the data acquisition is on or off, 
                'wave_sel' :                                  indicates the detection mode, 0 for full wave mode, 1 for feature mode, 
                'channel_control' :                       indicates on/off of the individual channels (in hex form), 
            }
    """

    print('readSciConfig: reading science config file \"' + basic.extractFilename(filename) + '\"')
    configDict = {}
    try:
        with open(filename, 'r') as loadfin:
            configDict = json.load(loadfin)
    except:
        raise Exception('readSciConfig: error reading config file \"' + basic.extractFilename(filename) + '\"')
    return configDict

def extractSciRawData(rawData, udpPackPos, maxUdpReadout = -1, lastUdpPos = 0, udpPackLen = specificParam.udpPackLen):

    """
    Function for extracting raw science data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from rundata file\n
    udpPackPos : array-like, 
        positions of UDP packages in raw data\n
    maxUdpReadout : int, optional
        maximum number of UDP packages to readout in a single readout run, for data readout in multiple readout runs in order to process large \
            data packs, if `maxUdpReadout` <= 0 then the data will be readout in a single run (default readout process)\n
    lastUdpPos : int, optional
        position of last UDP pack read, given as the data pack's index in udpPackPos\n
    udpPackLen : int, optional
        length of UDP packages in rundata, in bytes
    
    Returns
    ----------
    sciRawData : bytes
        raw science data read from raw data
    """

    udpPackagesID = [] #ID of UDP package, to get the count of UDP packages lost
    sciRawDataList = []
    #Set end positions for current readout run
    if maxUdpReadout <= 0:
        #Single readout run, ends at position of last UDP package
        posEnd = len(udpPackPos)
    else:
        #Multiple readout runs, reads 'maxUdpReadout' UDP packages in every run until the position of last UDP package
        posEnd = min(lastUdpPos + maxUdpReadout, len(udpPackPos))
    #Extract science data in every UDP package
    for ipos in range(lastUdpPos, posEnd):
        sciRawDataList.extend(rawData[udpPackPos[ipos] + 11:udpPackPos[ipos] + 11 + udpPackLen])
        udpPackagesID.append(rawData[udpPackPos[ipos] + 8] * 256 + rawData[udpPackPos[ipos] + 9])
        if (len(udpPackagesID) > 2 and (udpPackagesID[-1] - udpPackagesID[-2]) > 100):
            #Too many UDP packages lost, indicates error in the first data, therefore previously read data are discarded
            sciRawDataList.clear()
            udpPackagesID.clear()
    #Check UDP package loss
    udpPackageGap = np.uint16(udpPackagesID[1:]) - np.uint16(udpPackagesID[:-1])
    udpPackageLoss = sum(udpPackageGap) - len(udpPackageGap)
    if udpPackageLoss < 0:
        udpPackageLoss = 0
    sciRawData = bytes(sciRawDataList)
    sciRawDataList.clear()
    print('extractSciRawData: Lost {} / {} udp packages'.format(udpPackageLoss, len(udpPackageGap)))
    
    return sciRawData
''''
def extractUdp(rawData):
    patternU = rb'\x1a\xcf\xfc\x1d\x00\x00\x00\x00.{4}\x00\x00\x00'
    udpPattern = re.compile(patternU, re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)
    rawDataList = []
    for pos in udpPos:
        lenthD = struct.unpack('>l', rawData[pos+28:pos+32])[0]
        rawDataList.extend(rawData[pos+32:pos+32+lenthD])
        print(lenthD)
    rawData = bytes(rawDataList)
    print(rawData)
    rawDataList.clear()
    return rawData
'''
def extractHKData_03b(rawData, hkPackLen = 52):

    """
    Function to extract HK data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from HK data file\n
    hkPackLen : int, optional
        length of UDP packages of HK data, in bytes

    Returns
    ----------
    hkData : dict
        extracted science data, in the form of a dict as: \n
        hkData = {
            'iMon':                      monitored SiPM current,
            'bias':                       SiPM bias,
            'temp':                     SiPM temperature,
            'timestamp':             timestamp of current HK package, 
            'iSys':                       system current, 
        }
    """

    #Find positions of UDP packages
    udpPattern = re.compile(specificParam.pattrenRef['HK_03b'], re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)

    #HK data to extract
    hkData = {}
    bias = []
    iMon = []
    temp = []
    timestamp = []
    iSys = []
    for ich in range(4):
        bias.append([])
        iMon.append([])
        temp.append([])
        iSys.append([])

    #Extract data from all packages
    channelLookup = [0, 3, 2, 1]
    for i in np.arange(len(udpPos)):
        HKdata = rawData[udpPos[i] + 12:udpPos[i] + 12 + hkPackLen]
        for ich in np.arange(4):
            bias[channelLookup[ich]].append(struct.unpack('>H', HKdata[2 * ich + 0:2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('>H', HKdata[8 + 2 * ich:8 + 2 * ich + 2])[0])
            curtemp = struct.unpack('>H', HKdata[16 + 2 * ich:16 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('>H', HKdata[32 + 2 * ich:32 + 2 * ich + 2])[0])
        timestamp.append(struct.unpack('>Q', HKdata[24:24 + 8])[0] / specificParam.internalFreq)

    #Transform data into physical values
    hkData = {
        'iMon':                 np.array(iMon) / 2 ** 12 * 2.5 / (1 + 49.9 / 499) / 499 * 1E6, 
        'bias':                  np.array(bias) / 2 ** 12 * 2.5 / (51.1 / (1000 + 51.1)), 
        'temp':                np.array(temp) / 2 ** 4 * 0.0625, 
        'timestamp':        np.array(timestamp) * 100., 
        'iSys':                  np.array(iSys) / 2 ** 12 * 2.5 / (0.05 * 4.7E3 / 100), 
    }
    hkData['bias'] = hkData['bias']- hkData['iMon']*499*1E-6
    return hkData

def extractHKData_normal(rawData, hkPackLen = specificParam.hkDataLen):

    """
    Function to extract HK data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from HK data file\n
    hkPackLen : int, optional
        length of UDP packages of HK data, in bytes

    Returns
    ----------
    hkData : dict
        extracted science data, in the form of a dict as: \n
        hkData = {
            'iMon':                      monitored SiPM current,
            'bias':                       SiPM bias,
            'temp':                     SiPM temperature,
            'timestamp':             timestamp of current HK package, 
            'iSys':                       system current, 
        }
    """
    
    #rawData = extractUdp(rawData)
                
    #Find positions of UDP packages
    # udpPattern = re.compile(specificParam.pattrenRef['HK'], re.S)
    udpPattern = re.compile(specificParam.pattrenRef['HK_new'], re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)

    #HK data to extract
    hkData = {}
    bias = []
    iMon = []
    temp = []
    timestamp = []
    iSys = []
    for ich in range(4):
        bias.append([])
        iMon.append([])
        temp.append([])
        iSys.append([])

    #Extract data from all packages
    #channelLookup = [0, 3, 2, 1] #Special for 03 before 202105
    channelLookup = [0, 1, 2, 3]
    for i in np.arange(len(udpPos)):
        #HKdata = rawData[udpPos[i] + 12:udpPos[i] + 12 + hkPackLen]
        HKdata = rawData[udpPos[i] + 36:udpPos[i] + 36 + hkPackLen]
        for ich in np.arange(4):
            '''
            bias[channelLookup[ich]].append(struct.unpack('>H', HKdata[2 * ich + 0:2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('>H', HKdata[8 + 2 * ich:8 + 2 * ich + 2])[0])
            curtemp = struct.unpack('>H', HKdata[16 + 2 * ich:16 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('>H', HKdata[32 + 2 * ich:32 + 2 * ich + 2])[0])
            '''
            bias[channelLookup[ich]].append(struct.unpack('>H', HKdata[23 + 2 * ich:23 + 2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('>H', HKdata[31 + 2 * ich:31 + 2 * ich + 2])[0])
            curtemp = struct.unpack('>H', HKdata[47 + 2 * ich:47 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('>H', HKdata[63 + 2 * ich:63 + 2 * ich + 2])[0])
        #timestamp.append(struct.unpack('>Q', HKdata[24:24 + 8])[0] / specificParam.internalFreq)
        timestamp.append(struct.unpack('>Q', HKdata[55:55 + 8])[0] / specificParam.internalFreq)

    #Transform data into physical values
    hkData = {
        'iMon':                 np.array(iMon) / 2 ** 12 * 2.5 / (1 + 49.9 / 499) / 499 * 1E6, 
        'bias':                  np.array(bias) / 2 ** 12 * 2.5 / (51.1 / (1000 + 51.1)), 
        'temp':                np.array(temp) / 2 ** 4 * 0.0625, 
        'timestamp':        np.array(timestamp) * 100., 
        'iSys':                  np.array(iSys) / 2 ** 12 * 2.5 / (0.05 * 4.7E3 / 100), 
    }
    hkData['bias'] = hkData['bias']- hkData['iMon']*499*1E-6
    return hkData
def extractHKData(rawData, hkPackLen = specificParam.hkDataLen):

    """
    Function to extract HK data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from HK data file\n
    hkPackLen : int, optional
        length of UDP packages of HK data, in bytes

    Returns
    ----------
    hkData : dict
        extracted science data, in the form of a dict as: \n
        hkData = {
            'iMon':                      monitored SiPM current,
            'bias':                       SiPM bias,
            'temp':                     SiPM temperature,
            'timestamp':             timestamp of current HK package, 
            'iSys':                       system current, 
        }
    """
    
    #rawData = extractUdp(rawData)
                
    #Find positions of UDP packages
    #udpPattern = re.compile(specificParam.pattrenRef['HK'], re.S)
    udpPattern = re.compile(specificParam.pattrenRef['HK_new'], re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)

    #HK data to extract
    hkData = {}
    bias = []
    iMon = []
    temp = []
    timestamp = []
    iSys = []
    for ich in range(4):
        bias.append([])
        iMon.append([])
        temp.append([])
        iSys.append([])

    #Extract data from all packages
    #channelLookup = [0, 3, 2, 1] #Special for 03 before 202105
    channelLookup = [0, 1, 2, 3]
    for i in np.arange(len(udpPos)):
        #HKdata = rawData[udpPos[i] + 12:udpPos[i] + 12 + hkPackLen]
        HKdata = rawData[udpPos[i] + 36:udpPos[i] + 36 + hkPackLen]
        for ich in np.arange(4):
            '''
            bias[channelLookup[ich]].append(struct.unpack('>H', HKdata[2 * ich + 0:2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('>H', HKdata[8 + 2 * ich:8 + 2 * ich + 2])[0])
            curtemp = struct.unpack('>H', HKdata[16 + 2 * ich:16 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('>H', HKdata[32 + 2 * ich:32 + 2 * ich + 2])[0])
            '''

            '''
            bias[channelLookup[ich]].append(struct.unpack('>H', HKdata[23 + 2 * ich:23 + 2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('>H', HKdata[31 + 2 * ich:31 + 2 * ich + 2])[0])
            curtemp = struct.unpack('>H', HKdata[47 + 2 * ich:47 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('>H', HKdata[63 + 2 * ich:63 + 2 * ich + 2])[0])
            '''
            
            bias[channelLookup[ich]].append(struct.unpack('<H', HKdata[23 + 2 * ich:23 + 2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('<H', HKdata[31 + 2 * ich:31 + 2 * ich + 2])[0])
            curtemp = struct.unpack('<H', HKdata[47 + 2 * ich:47 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('<H', HKdata[63 + 2 * ich:63 + 2 * ich + 2])[0])

        #timestamp.append(struct.unpack('>Q', HKdata[24:24 + 8])[0] / specificParam.internalFreq)
        #timestamp.append(struct.unpack('>Q', HKdata[55:55 + 8])[0] / specificParam.internalFreq)
        timestamp.append(struct.unpack('<Q', HKdata[55:55 + 8])[0] / specificParam.internalFreq)

    #Transform data into physical values
    hkData = {
        'iMon':                 np.array(iMon) / 2 ** 12 * 2.5 / (1 + 49.9 / 499) / 499 * 1E6, 
        'bias':                  np.array(bias) / 2 ** 12 * 2.5 / (51.1 / (1000 + 51.1)), 
        'temp':                np.array(temp) / 2 ** 4 * 0.0625, 
        'timestamp':        np.array(timestamp) * 100., 
        'iSys':                  np.array(iSys) / 2 ** 12 * 2.5 / (0.05 * 4.7E3 / 100), 
    }
    hkData['bias'] = hkData['bias']- hkData['iMon']*499*1E-6
    return hkData

def extractTimelineData(rawData, tlPackLen = specificParam.timelineDataLen):

    """
    Function to extract timeline data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from timeline data file\n
    hkPackLen : int, optional
        length of UDP packages of timeline data, in bytes

    Returns
    ----------
    hkData : dict
        extracted science data, in the form of a dict as: \n
        hkData = {
            'utc':                        extracted UTC value,
            'pps':                        extracted PPS value,
            'timestamp':             timestamp of current timeline package, 
        }
    """

    #rawData = extractUdp(rawData)

    #Find positions of UDP packages
    udpPattern = re.compile(specificParam.pattrenRef['time_new'], re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)

    #Timeline data to extract
    tlData = {}
    utc = []
    pps = []
    timestamp = []

    #Extract data from all packages
    for i in np.arange(len(udpPos)):
        #TLdata = rawData[udpPos[i] + 11:udpPos[i] + 11 + tlPackLen]
        TLdata = rawData[udpPos[i] + 36:udpPos[i] + 36 + tlPackLen]
        utc.append(struct.unpack('>L', TLdata[:4])[0])
        pps.append(struct.unpack('>L', TLdata[4:8])[0])
        timestamp.append(struct.unpack('>Q', TLdata[8:16])[0] / specificParam.internalFreq)

    #Transform data into physical values
    tlData = {
        'utc':                   np.array(utc), 
        'pps':                   np.array(pps), 
        'timestamp':        np.array(timestamp) * 100., 
    }

    return tlData
def extractTimelineData_03b(rawData, tlPackLen = specificParam.timelineDataLen):

    """
    Function to extract timeline data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from timeline data file\n
    hkPackLen : int, optional
        length of UDP packages of timeline data, in bytes

    Returns
    ----------
    hkData : dict
        extracted science data, in the form of a dict as: \n
        hkData = {
            'utc':                        extracted UTC value,
            'pps':                        extracted PPS value,
            'timestamp':             timestamp of current timeline package, 
        }
    """

    #Find positions of UDP packages
    udpPattern = re.compile(specificParam.pattrenRef['timeline'], re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)

    #Timeline data to extract
    tlData = {}
    utc = []
    pps = []
    timestamp = []

    #Extract data from all packages
    for i in np.arange(len(udpPos)):
        TLdata = rawData[udpPos[i] + 11:udpPos[i] + 11 + tlPackLen]
        utc.append(struct.unpack('>L', TLdata[:4])[0])
        pps.append(struct.unpack('>L', TLdata[4:8])[0])
        timestamp.append(struct.unpack('>Q', TLdata[8:16])[0] / specificParam.internalFreq)

    #Transform data into physical values
    tlData = {
        'utc':                   np.array(utc), 
        'pps':                   np.array(pps), 
        'timestamp':        np.array(timestamp), 
    }

    return tlData
def extractIVData(rawData, IVPackLen = specificParam.IVPackLen):

    """
    Function to extract IV scan data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from HK data file\n
    IVPackLen : int, optional
        length of IV scan data, in bytes

    Returns
    ----------
    IVData : dict
        extracted IV scan data, in the form of a dict as: \n
        IVData = {
            'iMon':                      monitored SiPM current, 
            'vMon':                     monitored SiPM bias, 
        }
    """

    #Find positions of UDP packages
    udpPattern = re.compile(specificParam.pattrenRef['IV'], re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)

    #HK data to extract
    IVData = {}
    vMon = []
    iMon = []
    for ich in range(4):
        vMon.append([])
        iMon.append([])

    nPoints = int(IVPackLen / 24)
    #Extract data from all packages
    for i in np.arange(len(udpPos)):
        for ip in range(nPoints):
            IVdata = rawData[udpPos[i] + 11:udpPos[i] + 11 + IVPackLen]
            for ich in np.arange(4):
                iMon[ich].append(struct.unpack('>H', IVdata[24 * ip + 6 * ich:24 * ip + 6 * ich + 2])[0])
                vMon[ich].append(struct.unpack('<f', IVdata[24 * ip + 6 * ich + 2:24 * ip + 6 * ich + 6])[0])

    #Transform data into physical values
    IVData = {
        'iMon':                 np.array(iMon) / 2 ** 12 * 2.5 / (1 + 49.9 / 499) / 499 * 1E6, 
        'vMon':                np.array(vMon), 
    }

    return IVData

def extractVbrData(rawData, VbrPackLen = specificParam.VbrPackLen):

    """
    Function to extract Vbr scan data from raw data

    Parameters
    ----------
    rawData : bytes, 
        raw data read from HK data file\n
    VbrPackLen : int, optional
        length of Vbr scan data, in bytes

    Returns
    ----------
    VbrData : dict
        extracted Vbr scan data, in the form of a dict as: \n
        VbrData = {
            'iMon':                      monitored SiPM current, 
            'vMon':                     monitored SiPM bias, 
        }
    """

    #Find positions of UDP packages
    udpPattern = re.compile(specificParam.pattrenRef['Vbr'], re.S)
    udpPos = basic.findPackPos(rawData, udpPattern)

    #HK data to extract
    VbrData = {}
    vMon = []
    iMon = []
    for ich in range(4):
        vMon.append([])
        iMon.append([])

    nPoints = int(VbrPackLen / 24)
    #Extract data from all packages
    for i in np.arange(len(udpPos)):
        for ip in range(nPoints):
            VbrRawData = rawData[udpPos[i] + 11:udpPos[i] + 11 + VbrPackLen]
            for ich in np.arange(4):
                iMon[ich].append(struct.unpack('>H', VbrRawData[24 * ip + 6 * ich:24 * ip + 6 * ich + 2])[0])
                vMon[ich].append(struct.unpack('<f', VbrRawData[24 * ip + 6 * ich + 2:24 * ip + 6 * ich + 6])[0])

    #Transform data into physical values
    VbrData = {
        'iMon':                 np.array(iMon) / 2 ** 12 * 2.5 / (1 + 200E3 / 499) / 499 * 1E6, 
        'vMon':                np.array(vMon), 
    }

    return VbrData

def extractSciEvents(sciRawData, sampleLen = 256, baseLen = 1, newDatapack = False, featureMode = False, doDeduplicate = True, \
    featureEventNum = specificParam.featureEventNum, checkCrc = True):

    """
    Function to extract science data packs from raw science data

    Parameters
    ----------
    sciRawData : bytes
        raw science data read from raw data\n
    sampleLen : int, optional
        sample length for pulse mode data\n
    baseLen : int, optional
        length of baseline for pulse mode data\n
    newDatapack : boolean, optional
        indicates whether the data is produced with new datapack\n
    featureMode : boolean, optional
        indicates whether the input file is obtained with feature mode\n
    doDeduplicate : boolean, optional
        indicates whether the data will be deduplicated, for BINARY(.dat) data of new programme (6th ver.) only\n
    featureEventNum : int, optional
        number of events contained in one single feature mode data pack
    checkCrc : boolean, optional
        indicates whether CRC check will be performed for the data

    Returns
    ----------
    sciData : dict
        extracted science data, in the form of a dict as: \n
        sciData = {
            'timestampEvt':         event timestamp,
            'channel':                  channel of current event,
            'eventID':                  ID of current event,
            'charge':                   charge of current event (sum of pulse), for pulse mode only
            'amplitude':              amplitude of current event (max of pulse), for pulse mode only
            'meanBaseline':         mean baseline value (baseline at specific pulse position), for pulse mode only
            'baselineArray':         baseline data, for pulse mode only
            'pulseArray':              pulse data, for pulse mode only
        }
    """

    #Extract science packages
    nBytesData = 0
    if newDatapack:
        if featureMode:
            nBytesData = featureEventNum * 24 + 8
        else:
            nBytesData = sampleLen * 2 + 32
    else:
        nBytesData = 16 + sampleLen * 2 + 24
    if newDatapack:
        if featureMode:
            sciPattern = re.compile(specificParam.pattrenRef['sciFeatureNew'], re.S)
        else:
            sciPattern = re.compile(specificParam.pattrenRef['sciNewHead'] + bytes('{}'.format(nBytesData), encoding = 'utf-8') + \
                specificParam.pattrenRef['sciNewTail'], re.S)
    else:
        sciPattern = re.compile(specificParam.pattrenRef['sciHead'] + bytes('{}'.format(nBytesData), encoding = 'utf-8') + specificParam.pattrenRef['sciTail'], \
            re.S)
    sciPos, sciPackLen = basic.findPackPos(sciRawData, sciPattern, True)

    #Extract science data
    timestampEvt = []
    channel = []
    eventID = []
    charge = []
    amplitude = []
    meanBaseline = []
    baselineArray = []
    pulseArray = []
    crcErr = 0

    #Deduplicate
    pulsePackageAll = [sciRawData[sciPos[ipos]:sciPos[ipos] + sciPackLen[ipos]] for ipos in range(len(sciPos))]
    if doDeduplicate:
        pulsePackageAll = basic.deduplicate(pulsePackageAll)

    for ipack in range(len(pulsePackageAll)):
        try:
            pulsePackage = pulsePackageAll[ipack]
            # decode the event data package
            if newDatapack:
                curchannel = struct.unpack('>H', pulsePackage[8:10])[0]
                cursampleLen = struct.unpack('>H', pulsePackage[10:12])[0]
                cureventNumber = struct.unpack('>L', pulsePackage[12:16])[0]
                if featureMode:
                    crcCheck = True
                    curtimestampEvt = []
                    curdataSum = []
                    curdataMax = []
                    curdataBase = []
                    crcCal = pulsePackage[0:16]
                    for ievt in range(featureEventNum):
                        if checkCrc:
                            crc = struct.unpack('>H', pulsePackage[38 + ievt * 24:40 + ievt * 24])[0]
                            crcCal += pulsePackage[16 + ievt * 24:36 + ievt * 24]
                            if not basic.crcCheck(crcCal, pulsePackage[38 + ievt * 24:40 + ievt * 24]):
                                crcCheck = False
                                break
                        curtimestampEvt.append(struct.unpack('>Q', pulsePackage[16 + ievt * 24:24 + ievt * 24])[0] / specificParam.internalFreq)
                        curdataSum.append(struct.unpack('>Q', pulsePackage[24 + ievt * 24:32 + ievt * 24])[0])
                        curdataMax.append(struct.unpack('>H', pulsePackage[32 + ievt * 24:34 + ievt * 24])[0])
                        curdataBase.append(struct.unpack('>H', pulsePackage[34 + ievt * 24:36 + ievt * 24])[0])
                else:
                    curtimestampEvt = struct.unpack('>Q', pulsePackage[-32:-24])[0] / specificParam.internalFreq
                    curdataSum = struct.unpack('>Q', pulsePackage[-24:-16])[0]
                    curdataMax = struct.unpack('>H', pulsePackage[-16:-14])[0]
                    curdataBase = struct.unpack('>H', pulsePackage[-14:-12])[0]
                    crcCheck = basic.crcCheck(pulsePackage[0:-12], pulsePackage[-10:-8])
                    curdataDomain = list(pulsePackage[16:-32])
                    highBit = curdataDomain[0::2]
                    lowBit = curdataDomain[1::2]
                    curpulse = 256 * np.array(highBit) + np.array(lowBit)

                if crcCheck:
                    if featureMode:
                        for idata in range(len(curtimestampEvt)):
                            timestampEvt.append(curtimestampEvt[idata])
                            channel.append(curchannel)
                            eventID.append(cureventNumber)
                            charge.append(curdataSum[idata])
                            amplitude.append(curdataMax[idata])
                            meanBaseline.append(curdataBase[idata])
                    else:
                        timestampEvt.append(curtimestampEvt)
                        channel.append(curchannel)
                        eventID.append(cureventNumber)
                        charge.append(curdataSum)
                        amplitude.append(curdataMax)
                        meanBaseline.append(curdataBase)
                        baselineArray.append(curpulse[0:baseLen])
                        pulseArray.append(curpulse)
                else:
                    crcErr += 1

            else:
                curchannel = struct.unpack('>L', pulsePackage[8:12])[0]
                cureventNumber = struct.unpack('>L', pulsePackage[12:16])[0]
                curtimestampEvt = struct.unpack('>Q', pulsePackage[16:24])[0] / specificParam.internalFreq
                cursampleLen = struct.unpack('>H', pulsePackage[30:32])[0]
                curdataSum = struct.unpack('>Q', pulsePackage[-24:-16])[0]
                curdataMax = struct.unpack('>H', pulsePackage[-14:-12])[0]
                curdataBase = struct.unpack('>H', pulsePackage[-10:-8])[0]
                curdataDomain = list(pulsePackage[32:-24])
                highBit = curdataDomain[0::2]
                lowBit = curdataDomain[1::2]
                curpulse = 256 * np.array(highBit) + np.array(lowBit)

                timestampEvt.append(curtimestampEvt)
                channel.append(curchannel)
                eventID.append(cureventNumber)
                charge.append(curdataSum)
                amplitude.append(curdataMax)
                meanBaseline.append(curdataBase)
                baselineArray.append(curpulse[0:baseLen])
                pulseArray.append(curpulse)

        except Exception as e:
            print(e)
            pass
            
    if newDatapack:
        print('extractSciEvents: ' + str(crcErr) + ' of ' + str(len(sciPos)) + ' data packs with crc error')

    #Pack the data into dictionary
    if len(sciPos) > 0:
        sciData = {
            'timestampEvt':         np.array(timestampEvt, dtype = np.float64), 
            'channel':                  np.array(channel, dtype = np.uint8), 
            'eventID':                  np.array(eventID, dtype = np.uint32), 
            'charge':                   np.array(charge, dtype = np.uint64), 
            'amplitude':              np.array(amplitude, dtype = np.uint16), 
            'meanBaseline':         np.array(meanBaseline, dtype = np.uint16), 
            'baselineArray':         np.array(baselineArray, dtype = np.uint16), 
            'pulseArray':              np.array(pulseArray, dtype = np.uint16), 
        }
    else:
        sciData = {}

    return sciData

def dataReadout(filename, configFile = '', singleFile = False, featureMode = False, newDatapack = False, timeCut = None, doDeduplicate = False, \
    useFitBaseline = True, maxUdpReadout = -1, reqNum = [], ignoreCrc = False, plot = False, energyCut = None, baseLen = specificParam.baseLen, \
    noUdp = False,ending="x_ray"):

    """
    Function for reading out single GRID03 raw outout file, including: \n
    a) Seperate data files, with rundata, HK and timeline data in different binary (.dat) files\n
    b) Single data files, with rundata, HK and timeline data all in one single binary (.dat) file\n

    Parameters
    ----------
    filename : str,
        name of the output file, or name of the rundata file when data is stored in seperate files, with full path included\n
    configFile : str, optional
        name of the config file, if not given, science config will be read using default config file name\n
    singleFile : boolean, optional
        indicates whether the data is stored in seperate files (rundata, HK, timeline) or single file\n
    featureMode : boolean, optional
        indicates whether the input file is obtained with feature mode\n
    newDatapack : boolean, optional
        indicates whether the data is produced with new datapack\n
    timeCut : int or list, optional
        cut of time data in seconds\n
    doDeduplicate : boolean, optional
        indicates whether the data will be deduplicated\n
    useFitBaseline : boolean, optional
        indicates whether the amplitude data will be calculated with measured mean baseline, if False the baseline value in the config file will be used\n
    maxUdpReadout : int, optional
        maximum number of UDP packages to readout in a single readout run, for data readout in multiple readout runs in order to process large \
            data packs, if `maxUdpReadout` <= 0 then the data will be readout in a single run (default readout process)\n
    reqNum: list, optional
        required run numbers for science and telemetry data, in the form of [telreq, scireq]\n
    ignoreCrc : boolean, optional
        indicates whether CRC check will be ignored for current file\n
    plot : boolean, optional
        indicates whether figures of some inner-processing data (namely baseline) will be plotted\n
    energyCut : int or list, optional
        cut of energy in keV, specially designed for some in-orbit data\n
    baseLen : int, optional
        number of points in pulse array to use for baseline fit, e.g. the first N points where N is the value of `baseLen`\n
    noUdp : boolean, optional
        indicates whether the science raw data is arranged in UDP-level packages\n
    ending : str
        "x_ray" for extractHKData for_Xray;\n
        "normal" for extractHKData for normal data;\n
        "03b" for different channel loop;\n
    Returns
    ----------
    sciExtracted : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciExtracted = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
                    'eventID':                  array of arrays, extracted event ID, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 
                }
        All data are in the form of `array(array, array, array, array)` for single scan(including multiple runs)
    telExtracted : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telExtracted = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values, 
                    'iSys':                       array of arrays, extracted monitored system current, 
                    'timestamp':             array-like or array of arrays, extracted timestamp values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }
        Channel-seperated data are in the form of `array(array, array, array, array)` for single scan(including multiple runs)\n
        Shared data are in the form of `array` for single scan(including multiple runs)

    Raises
    ----------
    Exception
        when config file name is not given for non-seperate file data, or error reading data file(s)
    """

    #------------------------------------------------------------Basic settings------------------------------------------------------------
    filenameNoPath = basic.extractFilename(filename)
    if singleFile:
        if configFile == '':
            raise Exception('dataReadout: for non-seperate file data, name of config file should be specified')
    else:
        #Default config file name and HK data filename
        filePath = filename.split(filenameNoPath)[0]
        fileSplit = re.split('rundata|.dat', filenameNoPath)
        hkFile = filePath + fileSplit[0] + 'HK' + fileSplit[1] + '.dat'
        timelineFile = filePath + fileSplit[0] + 'TimeLine' + fileSplit[1] + '.dat'
        if configFile == '':
            if noUdp:
                configFile = filePath + fileSplit[0] + 'SciConfig' + fileSplit[1] + '.json'
            else:
                configFile = filePath + fileSplit[0] + 'scienceConfig' + fileSplit[1] + '.json'

    print('dataReadout: processing \"' + filenameNoPath + '\"')

    #Just in case that there is no science, HK or timeline data
    sciEmpty = False
    hkEmpty = False
    tlEmpty = False

    #Basic output
    amp = []
    dataMax = [] #dataMax = amplitude - baseline
    baseline = []
    timestampEvt = []
    baselineArray = []
    pulseArray = []
    eventID = []
    sciNum = []
    temp = []
    bias = []
    iMon = []
    timestamp = []
    iSys = []
    timestampTL = []
    utc = []
    telNum = []
    for ich in range(4):
        amp.append([])
        dataMax.append([])
        baseline.append([])
        timestampEvt.append([])
        baselineArray.append([])
        pulseArray.append([])
        eventID.append([])
        sciNum.append([])
        temp.append([])
        bias.append([])
        iMon.append([])
    
    section = 1 #telemetry run number(section)
    scisection = 1 #science run number(section)
    #Beginning and ends(in uscount values) of science and telemetry runs
    sciBegin = []
    sciEnd = []
    telBegin = []
    telEnd = []

    #Set required run for binary files
    telreq = -1
    scireq = -1
    if len(reqNum) == 2:
        telreq = reqNum[0]
        scireq = reqNum[1]
    elif filenameNoPath in specificParam.reqFileRef:
        telreq = specificParam.reqFileRef[filenameNoPath][0]
        scireq = specificParam.reqFileRef[filenameNoPath][1]
    
    #Read science config file
    configDict = readSciConfig(configFile, noUdp)

    #----------------------------------------------------------------Rundata readout----------------------------------------------------------------
    #For some files, ignore CRC check to ensure there is still telemetry data for further processing
    if filenameNoPath in specificParam.ignoreCrcRef:
        ignoreCrc = True
    
    reqRead = False #if required run specified, readout ends when required run is read
    rawData = bytes()
    sciRawData = bytes()
    #Read raw data from file
    try:
        with open(filename, 'rb') as fin:
            rawData = fin.read()
    except:
        raise Exception('dataReadout: error reading rundata from data file \"' + filenameNoPath +'\"')
    
    #Find UDP pack positions (pre-processing)
    udpPattern = re.compile(specificParam.pattrenRef['udp'], re.S)
    udpPackPos = basic.findPackPos(rawData, udpPattern)

    #Readout of science data packs
    lastPos = 0
    sampleLen = int(configDict['sampling_depth'] if noUdp else configDict['SamplingDepth'])
    if featureMode:
        sampleLen = 0
    totalReadoutRuns = int(float(len(udpPackPos)) / float(maxUdpReadout)) + 1 if maxUdpReadout >= 0 else 1
    try:
        for irun in range(totalReadoutRuns):
            print('Readout run ' + str(irun + 1) + ' of ' + str(totalReadoutRuns))
            if noUdp:
                sciRawData = copy(rawData)
            else:
                sciRawData = extractSciRawData(rawData, udpPackPos, maxUdpReadout, lastPos)
            sciData = extractSciEvents(sciRawData, sampleLen, baseLen, newDatapack, featureMode, doDeduplicate, checkCrc = not ignoreCrc)
            #Find decrease points in timestampEvt to split science runs
            if len(sciBegin) == 0:
                sciBegin.append(sciData['timestampEvt'][0])
            qSci = np.where(np.array(sciData['timestampEvt'])[:-1] > np.array(sciData['timestampEvt'])[1:] + specificParam.splitRunTime)[0]
            cursciNum = np.ones(len(sciData['timestampEvt'])) * scisection
            if len(qSci > 0):
                lastSciPos = 0
                for isci in range(len(qSci)):
                    scisection += 1
                    if scireq > -1 and scisection > scireq:
                        #Exit when the required run is read
                        reqRead = True
                    cursciNum[lastSciPos:qSci[isci] + 1] = scisection
                    lastSciPos = qSci[isci] + 1
                    sciEnd.append(sciData['timestampEvt'][qSci[isci]])
                    sciBegin.append(sciData['timestampEvt'][qSci[isci] + 1])
            emptyChannel = []
            #Find channels with no data
            for ich in range(4):
                if len(np.where(sciData['channel'] == ich)[0]) == 0:
                    emptyChannel.append(ich)
            for ich in range(4):
                if ich in emptyChannel:
                    continue
                if scireq < 0:
                    qCh = np.where(sciData['channel'] == ich)[0]
                else:
                    qCh = np.where(np.array(sciData['channel'] == ich) & np.array( cursciNum == scireq))[0]
                dataMax[ich].extend(list(sciData['amplitude'][qCh]))
                baseline[ich].extend(list(sciData['meanBaseline'][qCh]))
                timestampEvt[ich].extend(list(sciData['timestampEvt'][qCh]))
                eventID[ich].extend(list(sciData['eventID'][qCh]))
                sciNum[ich].extend(list(cursciNum[qCh]))
                if not featureMode:
                    baselineArray[ich].extend(list(sciData['baselineArray'][qCh]))
                    pulseArray[ich].extend(list(sciData['pulseArray'][qCh]))
            lastPos += maxUdpReadout
            if reqRead:
                break
        sciEnd.append(sciData['timestampEvt'][-1])
    except:
        print('dataReadout: no science data is read out. Please check the integrity or the pattern of the data pack')
        sciEmpty = True

    #Transform the data to ndarray(np.array)
    dataMax = np.array(dataMax, dtype = object)
    baseline = np.array(baseline, dtype = object)
    timestampEvt = np.array(timestampEvt, dtype = object)
    eventID = np.array(eventID, dtype = object)
    sciNum = np.array(sciNum, dtype = object)
    if not featureMode:
        baselineArray = np.array(baselineArray, dtype = object)
        pulseArray = np.array(pulseArray, dtype = object)

    #Fit baseline
    if not sciEmpty:
        baselineEdges = None
        if specificParam.useBaselineEdges:
            baselineEdges = np.arange(500, 2000, 5)
        if featureMode:
            offsetMeasured = fitBaseline(filenameNoPath, baseline, specEdges = baselineEdges, plot = plot)
        else:
            offsetMeasured = fitBaseline(filenameNoPath, baselineArray, specEdges = baselineEdges, plot = plot)

        #Calculate amplitude with config/measured baseline
        offset = None
        if noUdp:
            offset = baseline
        else:
            offset = [int(configDict['channel' + str(ich) + 'BaseLine']) for ich in range(4)]
        if useFitBaseline:
            offset = offsetMeasured
        for ich in range(4):
            if specificParam.useEventBaseline:
                amp[ich] = np.array(dataMax[ich]) - np.array(baseline[ich])
            elif newDatapack:
                amp[ich] = np.array(dataMax[ich]) - np.ones(len(baseline[ich])) * float(offset[ich])
            else:
                amp[ich] = np.array(dataMax[ich]) + np.array(baseline[ich]) - np.ones(len(baseline[ich])) * float(offset[ich])
        amp = np.array(amp, dtype = object)

    #Check pulse shape for pulse mode data
    if not featureMode and specificParam.checkPulseShape and not sciEmpty:
        plotPulseShape(filenameNoPath, pulseArray)

    #----------------------------------------------------------------HK data readout----------------------------------------------------------------

    #Read raw data from file
    hkRawData = bytes()
    if singleFile:
        hkRawData = rawData
    else:
        try:
            with open(hkFile, 'rb') as fin:
                hkRawData = fin.read()
        except:
            raise Exception('dataReadout: error reading HK data from data file \"' + filenameNoPath +'\"')

    hkExtracter = {
        'x_ray': extractHKData,
        'normal': extractHKData_normal,
        '03b': extractHKData_03b,
    }
    hkData = hkExtracter[ending](hkRawData)
    
    temp = hkData['temp']
    bias = hkData['bias']
    iMon = hkData['iMon']
    timestamp = hkData['timestamp']
    iSys = hkData['iSys']

    #Use try-except just in case there is no HK data
    try:
        #Find decrease points in timestampEvt to split telemetry runs
        telBegin.append(timestamp[0])
        qTel = np.where(np.array(timestamp)[:-1] > np.array(timestamp)[1:] + specificParam.splitRunTime)[0]
        telNum = np.ones(len(timestamp)) * section
        if len(qTel > 0):
            lastTelPos = 0
            for itel in range(len(qTel)):
                section += 1
                if telreq > -1 and section > telreq:
                    #Exit when the required run is read
                    break
                telNum[lastTelPos:qTel[itel] + 1] = section
                lastTelPos = qTel[itel] + 1
                telEnd.append(timestamp[qTel[itel]])
                telBegin.append(timestamp[qTel[itel] + 1])
        telEnd.append(timestamp[-1])

        #Select data according to required runs
        if telreq > -1:
            qTel = np.where(telNum == telreq)[0]
            timestamp = np.array(timestamp)[qTel]
            for ich in range(4):
                temp[ich] = np.array(temp[ich])[qTel]
                bias[ich] = np.array(bias[ich])[qTel]
                iMon[ich] = np.array(iMon[ich])[qTel]
                iSys[ich] = np.array(iSys[ich])[qTel]
    except:
        print('dataReadout: no HK data is read out. Please check the integrity or the pattern of the data pack')
        hkEmpty = True

    #Transform the data to ndarray(np.array)
    temp = np.array(temp)
    iMon = np.array(iMon)
    bias = np.array(bias)
    timestamp = np.array(timestamp)
    iSys = np.array(iSys)
    telNum = np.array(telNum)
    
    #-------------------------------------------------------------Timeline data readout-------------------------------------------------------------

    #Read raw data from file
    tlRawData = bytes()
    if singleFile:
        tlRawData = rawData
    else:
        try:
            with open(timelineFile, 'rb') as fin:
                tlRawData = fin.read()
        except:
            raise Exception('dataReadout: error reading HK data from data file \"' + filenameNoPath +'\"')
    if ending == '03b':
        tlData = extractTimelineData_03b(tlRawData)
    else:
        tlData = extractTimelineData(tlRawData)
    
    utc = tlData['utc']
    pps = tlData['pps']
    timestampTL = tlData['timestamp']

    try:
        #Calculate UTC values corresponding to the timestamp in HK data
        utc = getUTC(filenameNoPath, timestampTL, utc, timestamp, plot)

        #Transform the data to ndarray(np.array)
        utc = np.array(utc)
    except:
        print('dataReadout: no timeline data is read out. Please check the integrity or the pattern of the data pack')
        tlEmpty = True

    #----------------------------------------------------Readout data preliminary processing----------------------------------------------------

    #Display info of each telemetry and science run
    if not hkEmpty:
        print('Telemetry scan data: ')
        for il in range(len(telBegin)):
            print('Scan #' + str(il))
            print(telBegin[il], telEnd[il])

    if not sciEmpty:
        print('Science scan data: ')
        for il in range(len(sciBegin)):
            print('Scan #' + str(il))
            print(sciBegin[il], sciEnd[il])

    #Time cut
    if filenameNoPath in specificParam.cutFileRef:
        timeCut = specificParam.cutFileRef[filenameNoPath]
    if timeCut is not None:
        if not hkEmpty:
            q1 = None
            if isinstance(timeCut, list):
                # q1 = (timestamp >= timeCut[0]) * (timestamp <= timeCut[1])
                q1 = (timestamp >= 0) * (timestamp <= timeCut[1])
            else:
                q1 = (timestamp >= timeCut)
            timestamp = timestamp[q1]
            utc = utc[q1]
            temp = temp[:, q1]
            iMon = iMon[:, q1]
            bias = bias[:, q1]
            iSys = iSys[:, q1]
            telNum = telNum[q1]
        if not sciEmpty:
            for ich in range(4):
                q2 = None
                if isinstance(timeCut, list):
                    q2 = (np.array(timestampEvt[ich]) >= timeCut[0]) * (np.array(timestampEvt[ich]) <= timeCut[1])
                else:
                    q2 = np.array(timestampEvt[ich]) >= timeCut
                timestampEvt[ich] = np.array(timestampEvt[ich])[q2]
                amp[ich] = np.array(amp[ich])[q2]
                eventID[ich] = np.array(eventID[ich])[q2]
                sciNum[ich] = np.array(sciNum[ich])[q2]
    
    #Energy cut
    if energyCut is not None and not sciEmpty and not hkEmpty:
        curcorr = basic.tempBiasCorrection(temp, bias, ver = '03', style = parameters.tempBiasStyle)[0]
        curenergy = None
        for ich in range(4):
            curenergy = basic.ecCorrection(np.array(amp[ich]) * curcorr[ich], True, ich)
            if isinstance(energyCut, list):
                q1 = (np.array(curenergy) >= energyCut[0]) * (np.array(curenergy) <= energyCut[1])
            else:
                q1 = np.array(curenergy) >= energyCut
            timestampEvt[ich] = np.array(timestampEvt[ich])[q1]
            amp[ich] = np.array(amp[ich])[q1]
            eventID[ich] = np.array(eventID[ich])[q1]
            sciNum[ich] = np.array(sciNum[ich])[q1]

    #Check for error temperature data
    if parameters.tempCheck and not hkEmpty:
        qTemp = np.array(temp[0]) > -np.inf
        for ich in range(4):
            curtemp = list(temp[ich])
            curtemp.append(temp[ich][-1])
            curtemp = np.array(curtemp)
            qTemp *= (abs(curtemp[1:] - curtemp[:-1]) <= np.sqrt(np.std(temp[ich]) ** 2 + parameters.tempErrSys ** 2))
        timestamp = timestamp[qTemp]
        utc = utc[qTemp]
        temp = temp[:, qTemp]
        iMon = iMon[:, qTemp]
        bias = bias[:, qTemp]
        iSys = iSys[:, qTemp]
        telNum = telNum[qTemp]

    #Check for error UTC data
    if parameters.utcCheck and not tlEmpty:
        utcLaunch = basic.convertStrToTimestamp(specificParam.launchTime)
        qUtc = np.array(utc) > utcLaunch
        timestamp = timestamp[qUtc]
        utc = utc[qUtc]
        temp = temp[:, qUtc]
        iMon = iMon[:, qUtc]
        bias = bias[:, qUtc]
        iSys = iSys[:, qUtc]
        telNum = telNum[qUtc]

    #------------------------------------------------------------Data output------------------------------------------------------------
    print('Data readout of \"' + filenameNoPath + '\" complete')

    #Create output dictionaries
    sciExtracted = {
        'amp':                     amp, 
        'timestampEvt':       timestampEvt, 
        'eventID':                eventID, 
        'sciNum':                 sciNum, 
    }

    telExtracted = {
        'tempSipm':              temp, 
        'iMon':                      iMon, 
        'bias':                       bias, 
        'iSys':                       iSys, 
        'timestamp':             timestamp, 
        'utc':                        utc, 
        'telNum':                  telNum, 
    }

    #Output
    return sciExtracted, telExtracted

#****************************************************************************************************************************************************************
#***************************************************************Basic data process functions**************************************************************
#****************************************************************************************************************************************************************

def fitBaseline(filename, baselineArray, baselineRange = specificParam.baselineRange, binWidth = specificParam.baselineBinWidth, specEdges = None, \
    plot = False):

    """
    Function to plot and fit basiline distribution

    Parameters
    ----------
    filename : str, 
        name of the file\n
    baselineArray : array of arrays, 
        measured baseline values of all 4 channels, in the form of `[array, array, array, array]`\n
    baselineRange : list, optional
        fit and plot range of baseline data, in the form of `[lower, upper]` with 0 <= lower < upper\n
    binWidth : int, optional
        bin width when creating baseline histograms, being nonnegative integer\n
    specEdges : list or array-like or None, optional
        pre-defined spectrum edges, if not None, `specEdges` will be used instead of `binWidth` and `baselineRange`\n
    plot : boolean, optional
        indicates whether the fit baseline spectrum will be plotted\n
    
    Returns
    ----------
    offsetMeasured : list, 
        fit average baseline of all 4 channels, in the form of `[float, float, float, float]`
    
    Raises
    ----------
    Exception
        when baselineRange is not given in correct form or error saving figure
    """

    if specEdges is None and (len(baselineRange) != 2 or baselineRange[0] < 0 or baselineRange[1] <= baselineRange[0]):
        raise Exception('fitBaseline: fit range not given in correct form')
    
    if plot:
        fig = plt.figure(figsize = (12, 8))
        gs = gridspec.GridSpec(2, 2, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        if parameters.saveFigNoPlot:
            saveFigPath = ''
            filenameNoPath = basic.extractFilename(filename)
            filenameNoAppend = filenameNoPath[0:filenameNoPath.find(filenameNoPath.split('.')[-1]) - 1]
            saveFigPath = parameters.saveFigPath + '/'
            if parameters.saveFigOwnPath:
                saveFigPath += (filenameNoAppend + '/')
            if not os.path.isdir(saveFigPath):
                os.makedirs(saveFigPath)
    
    offsetMeasured = []
    
    baselineRangeCh = [copy(baselineRange), copy(baselineRange), copy(baselineRange), copy(baselineRange)]
    for ich in range(4):
        try:
            baselineMin = np.min(baselineArray[ich])
            baselineMax = np.max(baselineArray[ich])
            if baselineMin <= baselineRange[0] or baselineMax >= baselineRange[1]:
                if baselineMax - baselineMin < 1500:
                    baselineRangeCh[ich][0] = np.max([(baselineMin + baselineMax) / 2 - 750, 0])
                    baselineRangeCh[ich][1] = (baselineMin + baselineMax) / 2 + 750
                else:
                    baselineRangeCh[ich][0] = np.max([baselineMin - 200, 0])
                    baselineRangeCh[ich][1] = baselineMax + 200
        except:
            pass
    specBaseline, x  = basic.getSpectrum(baselineArray, binWidth = binWidth, specRange = baselineRangeCh, specEdges = specEdges)
    #Do plot and fit for all 4 channels
    for ich in range(4):
        try:
            if specEdges is None:
                y = specBaseline[ich] / binWidth
            else:
                binWidthEdges = specEdges[1:] - specEdges[:-1]
                y = specBaseline[ich] / binWidthEdges
            #Do gaussian fit
            out = basic.doFitGaussian(x[ich], y)
            center = out['fit_center']
            centerErr = out['fit_center_err']
            sigma = out['fit_sigma']
            sigmaErr = out['fit_sigma_err']
            amplitude = out['fit_amplitude']
            amplitudeErr = out['fit_amplitude_err']
            offsetMeasured.append(round(center))
            if plot:
                #Plot fit result
                ax = fig.add_subplot(gs[ich])
                yFit = basic.gaussianFunction([amplitude, center, sigma], x[ich])
                ax.step(x[ich], y, label = 'ch%d' % ich)
                ax.plot(x[ich], yFit, '-', label = 'fit')
                ax.legend(loc = 'best')
                ax.set_title('Baseline distribution,ch%d' % ich)
                ax.set_xlabel('ADC channel')
                ax.set_ylabel('counts')
                ax.set_ylim([0, 1.1 * max(y)])
                ax.tick_params(direction = 'out')
                ax.text(min(x[ich]) + 0.01 * (max(x[ich]) - min(x[ich])), max(y) - 0.1 * (max(y) - min(y)), '\nu = ${0:.2f}\pm {1:.2f}$\nFWHM = ${2:.2f}\pm {3:.2f}$\nA = '
                    '${4:.2f}\pm {5:.2f}$'.format(center, centerErr, 2.355 * sigma, 2.355 * sigmaErr, amplitude, amplitudeErr), rotation = 0, \
                    horizontalalignment = 'left', verticalalignment = 'top', multialignment = 'left', fontsize = 10)
        except:
            print('fitBaseline: channel %d has no fit baseline center' % ich)
            offsetMeasured.append(0.)
            pass
    
    #Save or plot figure
    if plot:
        if parameters.saveFigNoPlot:
            try:
                fig.savefig('{}/BaselineDistribution_{}.png'.format(saveFigPath, filenameNoAppend), dpi = 300,bbox_inches = 'tight')
                plt.close(fig)
            except:
                print('fitBaseline: error saving count rate figure for file \"' + filenameNoPath + '\"')
        else:
            plt.show()
    
    return offsetMeasured

def plotPulseShape(filename, pulseArray):

    """
    Function to plot pulse shape for pulse mode data

    Parameters
    ----------
    filename : str, 
        name of the file\n
    pulseArray : array of array, 
        sampled pulse data of all 4 channels\n

    Returns
    ----------
    Nothing

    Raises
    ----------
    Exception
        when error saving figure
    """
    
    saveFigPath = ''
    filenameNoPath = basic.extractFilename(filename)
    filenameNoAppend = filenameNoPath[0:filenameNoPath.find(filenameNoPath.split('.')[-1]) - 1]
    saveFigPath = parameters.saveFigPath + '/'
    if parameters.saveFigOwnPath:
        saveFigPath += (filenameNoAppend + '/')
    if not os.path.isdir(saveFigPath):
        os.makedirs(saveFigPath)

    fig = plt.figure(figsize = (12, 8))
    gs = gridspec.GridSpec(2, 2, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
    ipulse = specificParam.checkPulseNum
    for ich in range(4):
        try:
            ax = fig.add_subplot(gs[ich])
            curpulse = pulseArray[ich][ipulse]
            pulseTime = np.arange(len(curpulse))
            ax.step(pulseTime, curpulse, label = 'Pulse data')
            ax.set_xlabel('timestamp')
            ax.set_ylabel('amplitude')
            ax.set_title('ch{}'.format(ich))
        except:
            pass
    
    if parameters.saveFigNoPlot:
        try:
            fig.savefig('{}/PulseShape_{}.png'.format(saveFigPath, filenameNoAppend), dpi = 300, bbox_inches = 'tight')
            plt.close(fig)
        except:
            print('plotPulseShape: error saving count rate figure for file \"' + filenameNoPath + '\"')
    else:
        plt.show()

    return

def getUTC(filename, timestampTL, utc, timestamp, plot = False):

    """
    Function to calculate UTC values corresponding to timestamp in HK data

    Parameters
    ----------
    filename : str, 
        name of the file\n
    timestampTL : array-like, 
        timestamp in timeline data pack\n
    utc : array-like, 
        utc in timeline data pack, with the same shape as timestampTL\n
    timestamp : array-like, 
        timestamp in HK data pack\n
    plot: boolean, optional
        indicates whether the fit utc-timestampTL curve will be plotted\n

    Returns
    ----------
    utcFinal : array-like, 
        fit utc values corresponding to timestamp\n
    
    Raises
    ----------
    Exception
        when error saving figure
    """

    result = basic.doFitLin(timestampTL, utc)
    utcFinal = result['fit_a'] * np.array(timestamp) + result['fit_b']
    utcFinal = np.array(utcFinal, dtype = int)

    #Plot fit curve
    if plot:
        #Save figure settings
        saveFigPath = ''
        filenameNoPath = basic.extractFilename(filename)
        filenameNoAppend = filenameNoPath[0:filenameNoPath.find(filenameNoPath.split('.')[-1]) - 1]
        saveFigPath = parameters.saveFigPath + '/'
        if parameters.saveFigOwnPath:
            saveFigPath += (filenameNoAppend + '/')
        if not os.path.isdir(saveFigPath):
            os.makedirs(saveFigPath)
        
        timestampPlot = timestampTL
        fig = plt.figure(figsize = (12, 8))
        gs = gridspec.GridSpec(1, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        ax = fig.add_subplot(gs[0])
        ax.scatter(timestampPlot, utc, label = 'Raw data')
        ax.set_xlabel('timestamp/s')
        ax.set_ylabel('utc timestamp/s')
        ax.set_title('UTC-timestamp data')
        ax.grid()
        
        if parameters.saveFigNoPlot:
            try:
                fig.savefig('{}/utcFit_{}.png'.format(saveFigPath, filenameNoAppend), dpi = 300, bbox_inches = 'tight')
                plt.close(fig)
            except:
                print('getUTC: error saving count rate figure for file \"' + filenameNoPath + '\"')
        else:
            plt.show()

    return utcFinal

def splitData(filename, sciData, telData, singlech = False, channel = -1, splitByRun = True, timeCut = None, energyCut = None):

    """
    Function to split the data of multiple runs .dat files, according to given set bias in gridParameters or according to telemetry and science \
        run numbers, depending on the value of parameter `splitByRun`

    Parameters
    ----------
    filename : str,
        name of the file\n
    sciData : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciData = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
                    'eventID':                  array of arrays, extracted event ID, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 
                }
    telData : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telData = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values, 
                    'iSys':                       array of arrays, extracted monitored system current, 
                    'timestamp':             array-like or array of arrays, extracted timestamp values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }
    singlech : boolean, optional
        indicates whether the data selection will be done to only one single channel, for split by bias only\n
    channel : int, optional
        reference channel for single-channeled data selection, in range[0-3], for split by bias only\n
    splitByRun : boolean, optional
        indicates the way of splitting the data, True for splitting data according to given set bias in gridParameters, False for splitting data according to \
            telemetry and science run numbers
    timeCut : int or list, optional
        cut of time data in seconds, specially designed for in-orbit data\n
    energyCut : int or list, optional
        cut of energy in keV, specially designed for some in-orbit data\n

    Returns
    ----------
    sciDataFinal : dict
        split science data, in the form of a dictionary as\n
        sciDataFinal = {
                    'amp':                       array of arrays, split ADC amplitude, in the form of `array(array(array, len = NumOfRuns), \
                        array(array, len = NumOfRuns), array(array, len = NumOfRuns), array(array, len = NumOfRuns))`,
                    'timestampEvt':         array of arrays, split event uscount, in the same form as 'amp',
                    'eventID':                  array of arrays, split event ID, in the same form as 'amp', 
                    'sciNum':                   array of arrays, split science run number of each recorded event, 
                }
    telDataFinal : dict
        split telemetry data, in the form of a dictionary as\n
        telDataFinal = {
                    'tempSipm':              array of arrays, split SiPM temperature, in similar form as 'amp' in sciDataFinal, 
                    'iMon':                      array of arrays, split monitored current, in the same form as 'tempSipm', 
                    'bias':                       array of arrays, split bias, in the same form as 'tempSipm', 
                    'iSys':                       array of arrays, split monitored system current, in the same form as 'tempSipm', 
                    'timestamp':             array-like or array of arrays, split uscount values, in the form of `array(array, len = NumOfRuns)`, 
                    'utc':                        array-like or array of arrays, split UTC readout value, in the same form as 'uscount', 
                    'telNum':                  array-like or array of arrays, split telemetry run number of each temetry data, in the same form as 'timestamp', 
                }

    Raises
    ----------
    Execption
        when the data is single-channeled and given channel number is out of bound [0-3]
    """

    if singlech:
        if not basic.isChannel(channel):
            raise Exception('splitScans: channel number out of bound[0-3]')
    
    #Read each data from extracted data list
    amp = sciData['amp']
    timestampEvt = sciData['timestampEvt']
    eventID = sciData['eventID']
    sciNum = sciData['sciNum']
    tempSipm = telData['tempSipm']
    iMon = telData['iMon']
    bias = telData['bias']
    iSys = telData['iSys']
    timestamp = telData['timestamp']
    utc = telData['utc']
    telNum = telData['telNum']

    #Construct new data lists
    ampFinal = []
    timestampEvtFinal = []
    eventIDFinal = []
    sciNumFinal = []
    tempSipmFinal = []
    iMonFinal = []
    biasFinal = []
    iSysFinal = []
    timestampFinal = []
    utcFinal = []
    telNumFinal = []
    for ich in range(4):
        ampFinal.append([])
        timestampEvtFinal.append([])
        eventIDFinal.append([])
        sciNumFinal.append([])
        tempSipmFinal.append([])
        biasFinal.append([])
        iSysFinal.append([])
        iMonFinal.append([])
    
    #--------------------------------------------------Split data according to science and telemetry runs--------------------------------------------------
    if splitByRun:
        #Get starting and ending run numbers for the data cut
        curSciNum = -1
        curTelNum = -1
        endSciNum = -1
        endTelNum = -1
        print('splitScans: Splitting data according to science and telemetry runs')
        #Split telemetry run data
        if len(telNum) > 0:
            endTelNum = max(telNum)
            curTelNum = telNum[0]
            while curTelNum <= endTelNum:
                print('Processing telemetry scan, #{} of {}'.format(curTelNum, endTelNum))
                if filename in specificParam.unusedTelRef and curTelNum in specificParam.unusedTelRef[filename]:
                    curTelNum += 1
                    continue
                qTel = np.array(telNum) == curTelNum
                timestampFinal.append(list(np.array(timestamp)[qTel]))
                utcFinal.append(list(np.array(utc)[qTel]))
                telNumFinal.append(list(np.array(telNum)[qTel]))
                for ich in range(4):
                    tempSipmFinal[ich].append(list(np.array(tempSipm[ich])[qTel]))
                    iMonFinal[ich].append(list(np.array(iMon[ich])[qTel]))
                    iSysFinal[ich].append(list(np.array(iSys[ich])[qTel]))
                    biasFinal[ich].append(list(np.array(bias[ich])[qTel]))
                curTelNum += 1

        #Split science run data
        for ich in range(4):
            if len(sciNum[ich]) == 0:
                continue
            endSciNum = max(sciNum[ich])
            curSciNum = sciNum[ich][0]
            while curSciNum <= endSciNum:
                print('Processing science scan of channel{}, #{} of {}'.format(ich, curSciNum, endSciNum))
                if filename in specificParam.unusedSciRef and curSciNum in specificParam.unusedSciRef[filename]:
                    curSciNum += 1
                    continue
                qSci = np.array(sciNum[ich]) == curSciNum
                ampFinal[ich].append(list(np.array(amp[ich])[qSci]))
                timestampEvtFinal[ich].append(list(np.array(timestampEvt[ich])[qSci]))
                eventIDFinal[ich].append(list(np.array(eventID[ich])[qSci]))
                sciNumFinal[ich].append(list(np.array(sciNum[ich])[qSci]))
                curSciNum += 1
        
        #Time cut
        if timeCut is not None:
            #Telemetry data cut
            for itel in range(len(timestampFinal)):
                q1 = None
                if isinstance(timeCut, list):
                    q1 = (np.array(timestampFinal[itel]) >= timeCut[0]) * (np.array(timestampFinal[itel]) <= timeCut[1])
                else:
                    q1 = np.array(timestampFinal[itel]) >= timeCut
                    timestampFinal[itel] = list(np.array(timestampFinal[itel])[q1])
                    utcFinal[itel] = list(np.array(utcFinal[itel])[q1])
                    telNumFinal[itel] = list(np.array(telNumFinal[itel])[q1])
                    for ich in range(4):
                        tempSipmFinal[ich][itel] = list(np.array(tempSipmFinal[ich][itel])[q1])
                        iMonFinal[ich][itel] = list(np.array(iMonFinal[ich][itel])[q1])
                        iSysFinal[ich][itel] = list(np.array(iSysFinal[ich][itel])[q1])
                        biasFinal[ich][itel] = list(np.array(biasFinal[ich][itel])[q1])
            
            #Science data cut
            for ich in range(4):
                for isci in range(len(timestampEvtFinal[ich])):
                    q2 = None
                    if isinstance(timeCut, list):
                        q2 = (np.array(timestampEvtFinal[ich][isci]) >= timeCut[0]) * (np.array(timestampEvtFinal[ich][isci]) <= timeCut[1])
                    else:
                        q2 = np.array(timestampEvtFinal[ich, isci]) >= timeCut
                    ampFinal[ich][isci] = list(np.array(ampFinal[ich][isci])[q2])
                    timestampEvtFinal[ich][isci] = list(np.array(timestampEvtFinal[ich][isci])[q2])
                    eventIDFinal[ich][isci] = list(np.array(eventIDFinal[ich][isci])[q2])
                    sciNumFinal[ich][isci] = list(np.array(sciNumFinal[ich][isci])[q2])
    
        #Energy cut
        if energyCut is not None:
            tempCorr = np.array(tempSipmFinal, dtype = object)
            biasCorr = np.array(biasFinal, dtype = object)
            for ich in range(4):
                if len(timestampFinal) != len(timestampEvtFinal[ich]):
                    print('splitData: Unable to conduct energy cut due to non-corresponding number of science and telemetry runs. Cut session aborted')
                    break
                for isci in range(len(timestampEvtFinal[ich])):
                    curcorr = basic.tempBiasCorrection(tempCorr[:, isci], biasCorr[:, isci], ver = '03', style = parameters.tempBiasStyle)[0]
                    curenergy = None
                    for ich in range(4):
                        curenergy = basic.ecCorrection(np.array(ampFinal[ich][isci]) * curcorr[ich], True, ich)
                        if isinstance(energyCut, list):
                            q1 = (np.array(curenergy) >= energyCut[0]) * (np.array(curenergy) <= energyCut[1])
                        else:
                            q1 = np.array(curenergy) >= energyCut
                        timestampEvtFinal[ich][isci] = list(np.array(timestampEvtFinal[ich][isci])[q1])
                        ampFinal[ich][isci] = list(np.array(ampFinal[ich][isci])[q1])
                        eventIDFinal[ich][isci] = list(np.array(eventIDFinal[ich][isci])[q1])
                        sciNumFinal[ich][isci] = list(np.array(sciNumFinal[ich][isci])[q1])
    
    #--------------------------------------------------Split data according to specified bias values--------------------------------------------------
    else:
        #Get bias reference values and cut boundary for the data cut
        refVaules = specificParam.biasRef['bias_set']
        bound = specificParam.biasRef['bound']
        if filename in specificParam.unusedRef:
            #Skip unused values in some files
            refVaules = specificParam.unusedRef[filename]
        print('splitScans: Splitting data according to set bias ' + str(refVaules) + ' with bounds ' + str(bound))
        for biasSet in refVaules:
            timeBegin = -1
            timeEnd = np.max(timestamp) + 1
            #Record beginning and ending time of each split runs
            if singlech:
                qBias = (np.array(bias[channel]) >= (biasSet - bound)) * (np.array(bias[channel]) <= (biasSet + bound))
                if timeBegin < np.min(np.array(timestamp)[qBias]):
                    timeBegin = np.min(np.array(timestamp)[qBias])
                if timeEnd > np.max(np.array(timestamp)[qBias]):
                    timeEnd = np.max(np.array(timestamp)[qBias])
            else:
                for ich in range(4):
                    qBias = (np.array(bias[ich]) >= (biasSet - bound)) * (np.array(bias[ich]) <= (biasSet + bound))
                    if timeBegin < np.min(np.array(timestamp)[qBias]):
                        timeBegin = np.min(np.array(timestamp)[qBias])
                    if timeEnd > np.max(np.array(timestamp)[qBias]):
                        timeEnd = np.max(np.array(timestamp)[qBias])
            #Split the data according to the beginning and ending times
            qTel = (np.array(timestamp) >= timeBegin) * (np.array(timestamp) <= timeEnd)
            timestampFinal.append(list(np.array(timestamp)[qTel]))
            utcFinal.append(list(np.array(utc)[qTel]))
            telNumFinal.append(list(np.array(telNum)[qTel]))
            for ich in range(4):
                tempSipmFinal[ich].append(list(np.array(tempSipm[ich])[qTel]))
                iMonFinal[ich].append(list(np.array(iMon[ich])[qTel]))
                iSysFinal[ich].append(list(np.array(iSys[ich])[qTel]))
                biasFinal[ich].append(list(np.array(bias[ich])[qTel]))
                #Split science data
                qSci = (np.array(timestampEvt[ich]) >= timeBegin) * (np.array(timestampEvt[ich]) <= timeEnd)
                ampFinal[ich].append(np.array(amp[ich])[qSci])
                timestampEvtFinal[ich].append(list(np.array(timestampEvt[ich])[qSci]))
                eventIDFinal[ich].append(list(np.array(eventID[ich])[qSci]))
                sciNumFinal[ich].append(list(np.array(sciNum[ich])[qSci]))
    
    #Construct new data dictionaries with split data
    sciDataFinal = {
        'amp':                      np.array(ampFinal, dtype = object), 
        'timestampEvt':        np.array(timestampEvtFinal, dtype = object), 
        'eventID':                 np.array(eventIDFinal, dtype = object), 
        'sciNum':                  np.array(sciNumFinal, dtype = object), 
    }
    telDataFinal = {
        'tempSipm':              np.array(tempSipmFinal, dtype = object), 
        'iMon':                      np.array(iMonFinal, dtype = object), 
        'bias':                       np.array(biasFinal, dtype = object), 
        'iSys':                       np.array(iSysFinal, dtype = object), 
        'timestamp':             np.array(timestampFinal, dtype = object), 
        'utc':                        np.array(utcFinal, dtype = object), 
        'telNum':                  np.array(telNumFinal, dtype = object), 
    }
    return sciDataFinal, telDataFinal

def deleteEmptyScan(sciData, telData):

    """
    Auxiliary function to delete empty sublists for ranged multiple runs data

    Parameters
    ----------
    sciData : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciData = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
                    'eventID':                  array of arrays, extracted event ID, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 
                }
    telData : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telData = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values, 
                    'iSys':                       array of arrays, extracted monitored system current, 
                    'timestamp':             array-like or array of arrays, extracted timestamp values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }
    scanRange : list,
        list containing the scan range for multiple scans\n
    rateStyle : str,
        the style of calculating real count rate, '' for none, 's' for calculation with small data packs(512byte), 'l' for calculation with large data \
            packs(4096byte)\n
    newProgramme : boolean, optional
        indicates whether the data comes from new hardware programme(6th ver.)\n

    Returns
    ----------
    sciDataCleaned : dict
        science data with empty runs removed, in the form of a dictionary as\n
        sciDataCleaned = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
                    'eventID':                  array of arrays, extracted event ID, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 
                }
    telDataCleaned : dict
        telemetry data with empty runs removed, in the form of a dictionary as\n
        telDataCleaned = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values, 
                    'iSys':                       array of arrays, extracted monitored system current, 
                    'timestamp':             array-like or array of arrays, extracted timestamp values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }
    """

    #Read each data from extracted data list
    amp = sciData['amp']
    timestampEvt = sciData['timestampEvt']
    eventID = sciData['eventID']
    sciNum = sciData['sciNum']
    tempSipm = telData['tempSipm']
    iMon = telData['iMon']
    bias = telData['bias']
    iSys = telData['iSys']
    timestamp = telData['timestamp']
    utc = telData['utc']
    telNum = telData['telNum']

    nScan = len(timestamp)

    #Convert all data to list of lists in order to delete empty scans
    amp = list(amp)
    timestampEvt = list(timestampEvt)
    eventID = list(eventID)
    sciNum = list(sciNum)
    tempSipm = list(tempSipm)
    iMon = list(iMon)
    iSys = list(iSys)
    bias = list(bias)
    timestamp = list(timestamp)
    utc = list(utc)
    telNum = list(telNum)
    for isc in range(nScan):
        timestamp[isc] = list(timestamp[isc])
        utc[isc] = list(utc[isc])
        telNum[isc] = list(telNum[isc])
    for ich in range(4):
        amp[ich] = list(amp[ich])
        timestampEvt[ich] = list(timestampEvt[ich])
        eventID[ich] = list(eventID[ich])
        sciNum[ich] = list(sciNum[ich])
        tempSipm[ich] = list(tempSipm[ich])
        iMon[ich] = list(iMon[ich])
        iSys[ich] = list(iSys[ich])
        bias[ich] = list(bias[ich])
        for isc in range(nScan):
            amp[ich][isc] = list(amp[ich][isc])
            timestampEvt[ich][isc] = list(timestampEvt[ich][isc])
            eventID[ich][isc] = list(eventID[ich][isc])
            sciNum[ich][isc] = list(sciNum[ich][isc])
            tempSipm[ich][isc] = list(tempSipm[ich][isc])
            iMon[ich][isc] = list(iMon[ich][isc])
            iSys[ich][isc] = list(iSys[ich][isc])
            bias[ich][isc] = list(bias[ich][isc])

    #Delete empty runs
    currPos = 0
    while currPos < len(timestamp):
        if len(timestamp[currPos]) == 0:
            del timestamp[currPos], utc[currPos], telNum[currPos]
            for ich in range(4):
                del amp[ich][currPos], timestampEvt[ich][currPos], eventID[ich][currPos], sciNum[ich][currPos], tempSipm[ich][currPos], iMon[ich][currPos], \
                    iSys[ich][currPos], bias[ich][currPos]
        else:
            currPos += 1

    #Convert data back to ndarray form
    sciDataCleaned = {
        'amp':                      np.array(amp, dtype = object), 
        'timestampEvt':        np.array(timestampEvt, dtype = object), 
        'eventID':                 np.array(eventID, dtype = object), 
        'sciNum':                  np.array(sciNum, dtype = object), 
    }
    telDataCleaned = {
        'tempSipm':              np.array(tempSipm, dtype = object), 
        'iMon':                      np.array(iMon, dtype = object), 
        'iSys':                       np.array(iSys, dtype = object), 
        'bias':                       np.array(bias, dtype = object), 
        'timestamp':             np.array(timestamp, dtype = object), 
        'utc':                        np.array(utc, dtype = object), 
        'telNum':                  np.array(telNum, dtype = object), 
    }

    #Thanks to ndarray, ALL THESE matter is needed to delete a SINGLE sublist!
    return sciDataCleaned, telDataCleaned

def addDictToExistingData(newSciData, sciData, newTelData, telData, isCi = 0):

    """
    Auxiliary function to add new readout data (in the form of dict) to existing readout data

    Parameters
    ----------
    newSciData : dict
        dictionary of new readout science data\n
    sciData : dict
        dictionary of existing science data, possibly missing essential keys\n
    newTelData : dict
        dictionary of new readout telemetry data\n
    telData : dict
        dictionary of existing telemetry data, possibly missing essential keys\n
    isCi : int, optional
        indicates whether the input file has multiple scans/runs, 0 and 1 for single scan/run, 2 for multiple scans/runs\n
    
    Returns
    ----------
    sciData : dict
        dictionary of existing science data with new data added\n
    telData : dict
        dictionary of existing telemetry data with new data added\n
    """

    #Process science data
    sciMultiChannel = ['amp', 'timestampEvt', 'eventID', 'sciNum']
    newData = False
    for key in newSciData:
        if key not in sciData:
            sciData.update({
                key:        newSciData[key], 
            })
            if isCi == 2:
                sciData.update({
                    'multiScan':        True, 
                })
        else:
            newData = True
            if key in sciMultiChannel:
                if not sciData['multiScan']:
                    newKeyData = []
                    for ich in range(4):
                        newKeyData.append([sciData[key][ich]])
                    sciData[key] = newKeyData
                sciData[key] = list(sciData[key])
                for ich in range(4):
                    sciData[key][ich] = list(sciData[key][ich])
                    if isCi == 2:
                        for isc in range(len(newSciData[key][ich])):
                            sciData[key][ich].append(newSciData[key][ich][isc])
                    else:
                        sciData[key][ich].append(newSciData[key][ich])
                sciData[key] = np.array(sciData[key], dtype = object)
            else:
                if not sciData['multiScan']:
                    sciData[key] = [sciData[key]]
                sciData[key] = list(sciData[key])
                if isCi == 2:
                    for isc in range(len(newSciData[key])):
                        sciData[key].append(newSciData[key][isc])
                else:
                    sciData[key].append(newSciData[key])
                sciData[key] = np.array(sciData[key], dtype = object)
    if newData:
        sciData.update({
            'multiScan':        True, 
        })
    
    #Process telemetry data
    telMultiChannel = ['tempSipm', 'iMon', 'bias', 'iSys']
    newData = False
    for key in newTelData:
        if key not in telData:
            telData.update({
                key:        newTelData[key], 
            })
            if isCi == 2:
                telData.update({
                    'multiScan':        True, 
                })
        else:
            newData = True
            if key in telMultiChannel:
                if not telData['multiScan']:
                    newKeyData = []
                    for ich in range(4):
                        newKeyData.append([telData[key][ich]])
                    telData[key] = newKeyData
                telData[key] = list(telData[key])
                for ich in range(4):
                    telData[key][ich] = list(telData[key][ich])
                    if isCi == 2:
                        for isc in range(len(newTelData[key][ich])):
                            telData[key][ich].append(newTelData[key][ich][isc])
                    else:
                        telData[key][ich].append(newTelData[key][ich])
                telData[key] = np.array(telData[key], dtype = object)
            else:
                if not telData['multiScan']:
                    telData[key] = [telData[key]]
                telData[key] = list(telData[key])
                if isCi == 2:
                    for isc in range(len(newTelData[key])):
                        telData[key].append(newTelData[key][isc])
                else:
                    telData[key].append(newTelData[key])
                telData[key] = np.array(telData[key], dtype = object)
    if newData:
        telData.update({
            'multiScan':        True, 
        })
    
    #Output
    return sciData, telData

def fitRateCorrect(filename, timestampEvt, eventID, singlech = False, channel = -1, plot = True, odr = False):
    
    """
    Function for calculating correct total count rate and its error with the timeCorrect calculated when reading the data

    filename : str,
        name of the input data file\n
    timestampEvt : array of arrays,
        event timestamp extracted from science data\n
    eventID : array of arrays,
        event ID extracted from science data\n
    singlech : boolean, optional
        indicates whether the data selection will be done to only one single channel, for split by bias only\n
    channel : int, optional
        reference channel for single-channeled data selection, in range[0-3], for split by bias only\n
    plot : boolean, optional
        indicates the whether the figure of the fit result will be plotted\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n

    Returns
    ----------
    rateAll : list,
        correct total count rate of all 4 channels\n
    rateAllErr : list,
        error of correct total count rate of all 4 channels\n

    Raises
    ----------
    Exception
        when error saving figure
    """

    rateAll = []
    rateAllErr = []
    if not singlech:
        rateAll = []
        rateAllErr = []
    if plot:
        fig = plt.figure(figsize=(12, 8))
        if singlech:
            gs = gridspec.GridSpec(1, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        else:
            gs = gridspec.GridSpec(2, 2, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        
    for ich in range(4):
        if singlech and ich != channel:
            rateAll.append(1)
            rateAllErr.append(1)
            continue
        #Fit count rate for each channel
        timeData = (timestampEvt[ich] - timestampEvt[ich][0])
        countData = eventID[ich] - eventID[ich][0]
        result = basic.doFitLin(timeData, countData, odr = odr)
        a = result['fit_a']
        aErr = result['fit_a_err']
        b = result['fit_b']
        rateAll.append(a)
        rateAllErr.append(aErr)
        #Plot part
        if plot:
            fitCount = basic.linearFunction([a, b], timeData)
            ax = fig.add_subplot(gs[ich])
            ax.scatter(timeData, countData, label = 'raw data', zorder = 1, c = 'b')
            ax.plot(timeData, fitCount, label = 'Linear Fit', c = 'r')
            ax.set_ylim([0., 1.2 * max(np.max(countData), np.max(fitCount))])
            ax.text(np.average(timeData), max(fitCount) * 1.1, 'count rate = ' + str('%.2e' % a) + ' $\pm$ ' + str('%.3e' % aErr) + 'cps', fontsize = 10, \
                bbox = dict(facecolor = 'pink', alpha = 0.1), horizontalalignment = 'center', verticalalignment = 'center')
            ax.set_xlabel('timestamp/s')
            ax.set_ylabel('count')
            ax.legend(loc = 0)
            ax.grid()
    
    #Save/plot figure
    if parameters.saveFigNoPlot:
        filenameNoPath = basic.extractFilename(filename)
        filenameNoAppend = filenameNoPath[0:filenameNoPath.find(filenameNoPath.split('.')[-1]) - 1]
        savePath = parameters.saveFigPath + '/'
        if parameters.saveFigOwnPath:
            savePath += (filenameNoAppend + '/')
        if not os.path.isdir(savePath):
            os.makedirs(savePath)
        try:
            fig.savefig('{}/{}_rate_fit.png'.format(savePath, filenameNoAppend), transparent = True)
            plt.close(fig)
        except:
            print('fitRateCorrect: error saving count rate figure for file \"' + filenameNoPath + '\"')
    else:
        plt.show()
    
    return rateAll, rateAllErr

#****************************************************************************************************************************************************************
#***************************************************************Processed data I/O functions**************************************************************
#****************************************************************************************************************************************************************

def fileOutput(file, sciData, telData):

    """
    Function for writing designated readout data to files(numpy form binary file .npy)

    Parameters
    ----------
    file : str,
        name of the raw data file\n
    sciData : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciData = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
                    'eventID':                  array of arrays, extracted event ID, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 
                }
    telData : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telDataCleaned = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values, 
                    'iSys':                       array of arrays, extracted monitored system current, 
                    'timestamp':             array-like or array of arrays, extracted timestamp values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }

    Returns
    ----------
    boolean,
        True if all output files have been successfully written, False if not\n

    Notes
    ----------
    Output file naming rules(following the previous naming rules):\n
    (`channel` : channel number; `scancount` : scan number; `filename` : filename of raw data with affix[path] and suffix[\'.dat\'/\'.txt\'] removed)\n
    ----
    ADC amplitude: `amp_filename.npy`\n
    Event timestamp: `timestampEvt_filename.npy`\n
    Event ID: `eventID_filename.npy`\n
    Science run number: `sciNum_filename.npy`\n
    SiPM temperature: `tempSiPM_filename.npy`\n
    Monitored current:  `iMon_filename.npy`\n
    Bias: `bias_filename.npy`\n
    System current: `iSys_filename.npy`\n
    Timestamp: `timestamp_filename.npy`\n
    UTC: `utc_filename.npy`\n
    Telemetry run number: `telNum_filename.npy`\n
    """

    amp = sciData['amp']
    timestampEvt = sciData['timestampEvt']
    eventID = sciData['eventID']
    sciNum = sciData['sciNum']
    tempSipm = telData['tempSipm']
    iMon = telData['iMon']
    bias = telData['bias']
    iSys = telData['iSys']
    timestamp = telData['timestamp']
    utc = telData['utc']
    telNum = telData['telNum']

    print('fileOutput: writing output files of ' + file)
    outputSuccess = True
    errCount = 0

    #Read config file
    config = {}
    try:
        with open(parameters.outputConfig['03'], 'r') as fin:
            config = json.load(fin)
    except:
        print('fileOutput: unable to read file output configuration from file \"' + parameters.outputConfig['03'] + '\"')
        return False

    try:
        filename = file.split('.')[0]
        savePath = config['output_dir']
        if config['use_own_dir']:
            savePath += (filename + '/')
        if not os.path.exists(savePath):
            os.makedirs(savePath)

        #ADC amplitude
        if config['amp']:
            try:
                with open(savePath + 'amp_'+ filename + '.npy', 'wb') as foutamp:
                    np.save(foutamp, amp)
            except:
                print('fileOutput: Error writing ADC amplitude file \"' + 'amp_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Science run number
        if config['sciNum']:
            try:
                with open(savePath + 'sciNum_'+ filename + '.npy', 'wb') as foutsci:
                    np.save(foutsci, sciNum)
            except:
                print('fileOutput: Error writing science run number file \"' + 'sciNum_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #SiPM temperature
        if config['tempSipm']:
            try:
                with open(savePath + 'tempSiPM_' + filename + '.npy', 'wb') as fouttemp:
                    np.save(fouttemp, tempSipm)
            except:
                print('fileOutput: Error writing SiPM temperature file \"' + 'tempSiPM_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Monitored current
        if config['iMon']:
            try:
                with open(savePath + 'iMon_' + filename + '.npy', 'wb') as foutiMon:
                    np.save(foutiMon, iMon)
            except:
                print('fileOutput: Error writing monitored current file \"' + 'iMon_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Bias
        if config['bias']:
            try:
                with open(savePath + 'bias_' + filename + '.npy', 'wb') as foutbias:
                    np.save(foutbias, bias)
            except:
                print('fileOutput: Error writing bias file \"' + 'bias_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #System current
        if config['iSys']:
            try:
                with open(savePath + 'iSys_' + filename + '.npy', 'wb') as foutiSys:
                    np.save(foutiSys, iSys)
            except:
                print('fileOutput: Error writing system current file \"' + 'iSys_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Uscount
        if config['timestamp']:
            try:
                with open(savePath + 'timestamp_' + filename + '.npy', 'wb') as fouttime:
                    np.save(fouttime, timestamp)
            except:
                print('fileOutput: Error writing timestamp file \"' + 'timestamp_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #UTC
        if config['utc']:
            try:
                with open(savePath + 'utc_' + filename + '.npy', 'wb') as fouttime:
                    np.save(fouttime, utc)
            except:
                print('fileOutput: Error writing uscount file \"' + 'utc_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Telemetry run number
        if config['telNum']:
            try:
                with open(savePath + 'telNum_' + filename + '.npy', 'wb') as fouttel:
                    np.save(fouttel, telNum)
            except:
                print('fileOutput: Error writing telemetry run number file \"' + 'telNum_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Event uscount
        if config['timestampEvt']:
            try:
                with open(savePath + 'timestampEvt_' + filename + '.npy', 'wb') as fouttimeevt:
                    np.save(fouttimeevt, timestampEvt)
            except:
                print('fileOutput: Error writing event timestamp file \"' + 'timestampEvt_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Event ID
        if config['eventID']:
            try:
                with open(savePath + 'eventID_' + filename + '.npy', 'wb') as foutevtID:
                    np.save(foutevtID, eventID)
            except:
                print('fileOutput: Error writing event ID file \"' + 'eventID_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        print('File output complete with ' + str(errCount) + ' error(s)')

    except:
        print('fileOutput: error occurred while reading config file, abandoning file output process')
        outputSuccess = False

    return outputSuccess

def importData(filename, importPath, isCi = 0):

    """
    Function for importing data from output files

    Parameters
    ----------
    filename : list,
        names of the raw data with affix[path] and suffix[\'.dat\'/\'.txt\'] removed, currently supporting only numpy data(.npy) files\n
    importPath : list,
        paths of the import directories\n
    isCi : int, optional
        indicates whether the input file has CI part, with 0 for no CI, 1 for with CI\n

    Returns
    ----------
    sciExtracted : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciExtracted = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event timestamp values, which corresponds to all events recorded in amp,
                    'eventID':                  array of arrays, extracted event ID, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 
                }
        All data are in the form of `array(array, array, array, array)` for single scan(including multiple runs), or `array(array(array, len = NumOfScans), \
            array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans))` for multiple scans
    telExtracted : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telExtracted = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values, 
                    'iSys':                       array of arrays, extracted monitored system current, 
                    'timestamp':             array-like or array of arrays, extracted timestamp values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }
        Channel-seperated data are in the form of `array(array, array, array, array)` for single scan(including multiple runs), or \
            `array(array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans))` \
            for multiple scans\n
        Shared data are in the form of `array` for single scan(including multiple runs), or \
            `array(array, len = NumOfScans)` for multiple scans\n
    """

    #Read config file
    config = {}
    readConfig = True
    try:
        with open(specificParam.importConfig, 'r') as fin:
            config = json.load(fin)
    except:
        print('importData: unable to read file output configuration from file \"' + specificParam.importConfig + '\"')
        print('WARNING: empty data sets will be returned from \'importData\'')
        readConfig = False

    importNames = []
    #Preprocessing of file names and paths
    if len(importPath) == 0:
        importPath.append(config['import_dir'])
    if (len(filename) == 0 or len(filename) == 1 and filename[0] == '') and config['use_own_dir']:
        print('importData: unable to import from unspecified directory with \'use_own_dir\' == True in config file \"' + specificParam.importConfig + '\" and '
            'no file names given')
        print('WARNING: empty data sets will be returned from \'importData\'')
        readConfig = False
    elif len(filename) > 0:
        if len(filename) == 1 and filename[0] == '':
            filename = []
        else:
            for file in filename:
                rootname = file.spilt('/')[-1]
                importNames.append(rootname.split('.')[0])

    amp = [] #after CI for data with CI
    timestampEvt = []
    eventID = []
    sciNum = []
    tempSipm = []
    iMon = []
    bias = []
    iSys = []
    timestamp = [] #timestamp in telemetry data
    utc = []
    telNum = []

    #Create output dictionaries
    sciExtracted = {}
    telExtracted = {}

    for ich in range(4):
        amp.append([])
        timestampEvt.append([])
        eventID.append([])
        sciNum.append([])
        tempSipm.append([])
        iMon.append([])
        bias.append([])
        iSys.append([])

    errCount = 0
    if readConfig:
        for path in importPath:
            print('importData: reading files in ' + path)
            files = os.listdir(path)
            for file in files:
                #Select files to import
                extractFile = False
                if len(importNames) == 0:
                    #If no import name is specified, all numpy data(.npy) files within the directory are imported
                    if file.endswith('.npy'):
                        extractFile = True
                else:
                    #If import name(s) are specified, import only files containing import names
                    for name in importNames:
                        if file.endswith(name + '.npy'):
                            extractFile = True
                            break

                if extractFile:
                    #ADC amplitude files
                    if file.startswith('amp') and config['amp']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curamp = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curamp[ich])):
                                            amp[ich].append(curamp[ich][isc])
                                    else:
                                        amp[ich].append(curamp[ich])
                        except:
                            print('importData: error importing file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Science run number files
                    if file.startswith('sciNum') and config['sciNum']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                cursciNum = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(cursciNum[ich])):
                                            sciNum[ich].append(cursciNum[ich][isc])
                                    else:
                                        sciNum[ich].append(cursciNum[ich])
                        except:
                            print('importData: error importing file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Event uscount files
                    elif file.startswith('uscevt') and config['uscountEvt']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curtimestampEvt = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curtimestampEvt[ich])):
                                            timestampEvt[ich].append(curtimestampEvt[ich][isc])
                                    else:
                                        timestampEvt[ich].append(curtimestampEvt[ich])
                        except:
                            print('importData: error importing file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Sipm temperature files
                    elif file.startswith('tempSiPM') and config['tempSipm']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curtempSipm = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curtempSipm[ich])):
                                            tempSipm[ich].append(curtempSipm[ich][isc])
                                    else:
                                        tempSipm[ich].append(curtempSipm[ich])
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Monitored current files
                    elif file.startswith('iMon') and config['iMon']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curiMon = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curiMon[ich])):
                                            iMon[ich].append(curiMon[ich][isc])
                                    else:
                                        iMon[ich].append(curiMon[ich])
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #SiPM bias files
                    elif file.startswith('bias') and config['bias']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curbias = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curbias[ich])):
                                            bias[ich].append(curbias[ich][isc])
                                    else:
                                        bias[ich].append(curbias[ich])
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #System current files
                    elif file.startswith('iSys') and config['iSys']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curiSys = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curiSys[ich])):
                                            iSys[ich].append(curiSys[ich][isc])
                                    else:
                                        iSys[ich].append(curiSys[ich])
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Timestamp files
                    elif file.startswith('timestamp') and config['timestamp']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curtimestamp = np.load(fin, allow_pickle = True)
                                if isCi == 2:
                                    for isc in range(len(curtimestamp)):
                                        timestamp.append(curtimestamp[isc])
                                else:
                                    timestamp.append(curtimestamp)
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #UTC files
                    elif file.startswith('utc') and config['utc']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curutc = np.load(fin, allow_pickle = True)
                                if isCi == 2:
                                    for isc in range(len(curutc)):
                                        utc.append(curutc[isc])
                                else:
                                    utc.append(curutc)
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Telemetry run number files
                    elif file.startswith('telNum') and config['telNum']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curtelNum = np.load(fin, allow_pickle = True)
                                if isCi == 2:
                                    for isc in range(len(curtelNum)):
                                        telNum.append(curtelNum[isc])
                                else:
                                    telNum.append(curtelNum)
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

        print('Data import complete with ' + str(errCount) + ' error(s)')

    #Convert single scan data to correct form and create output dictionaries
    if config['timestamp']:
        if len(timestamp) == 1:
            timestamp = timestamp[0]
        telExtracted.update({
            'timestamp':              np.array(timestamp), 
        })
    if config['utc']:
        if len(utc) == 1:
            utc = utc[0]
        telExtracted.update({
            'utc':              np.array(utc), 
        })
    if config['amp']:
        for ich in range(4):
            if len(amp[ich]) == 1:
                amp[ich] = amp[ich][0]
        sciExtracted.update({
            'amp':              np.array(amp, dtype = object), 
        })
    if config['timestampEvt']:
        for ich in range(4):
            if len(timestampEvt[ich]) == 1:
                timestampEvt[ich] = timestampEvt[ich][0]
        sciExtracted.update({
            'timestampEvt':              np.array(timestampEvt, dtype = object), 
        })
    if config['tempSipm']:
        for ich in range(4):
            if len(tempSipm[ich]) == 1:
                tempSipm[ich] = tempSipm[ich][0]
        telExtracted.update({
            'tempSipm':              np.array(tempSipm), 
        })
    if config['iMon']:
        for ich in range(4):
            if len(iMon[ich]) == 1:
                iMon[ich] = iMon[ich][0]
        telExtracted.update({
            'iMon':              np.array(iMon), 
        })
    if config['bias']:
        for ich in range(4):
            if len(bias[ich]) == 1:
                bias[ich] = bias[ich][0]
        telExtracted.update({
            'bias':              np.array(bias), 
        })
    if config['iSys']:
        for ich in range(4):
            if len(iSys[ich]) == 1:
                iSys[ich] = iSys[ich][0]
        telExtracted.update({
            'iSys':              np.array(iSys), 
        })
    
    #Data output
    return sciExtracted, telExtracted
