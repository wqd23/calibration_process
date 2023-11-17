# -*- coding:utf-8 -*-

"""
GRID Basic Functions Library
----------

Basic functions for GRID data processing and fit\n
`v1.0.0` for GRID02 calibration result analysis, also quick view of GRID02 in-orbit data, by ydx and ghz\n
"""

#**************************************************************************************************************************************************
#***************************************************************Import packages***************************************************************
#**************************************************************************************************************************************************

from . import gridParametersCommon as parameters
import numpy as np
import lmfit
from scipy.odr import ODR, Model, RealData
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import struct
import crc16
import os
import re
import json, csv
import datetime
from copy import copy

#****************************************************************************************************************************************************************
#**********************************************************Basic data process and fit functions**********************************************************
#****************************************************************************************************************************************************************

#*******************************************************************************************************************************************************
#***************************************************************Auxiliariry functions***************************************************************
#*******************************************************************************************************************************************************

def isChannel(input):

    """
    Auxiliary function to determine whether the input is a channel number

    Parameters
    ----------
    input : str,
        the input string\n

    Returns
    ----------
    boolean,
        True if the input is a integer within [0-3], False if not
    """

    try:
        ch = int(input)
        if ch < 0 or ch > 3:
            raise Exception
    except:
        return False
    return True

def getSpectrum(amp, nbins = parameters.adcMax, singlech = False, specRange = [0., parameters.adcMax], specEdges = None, binWidth = None):

    """
    Function for getting the spectrum of the amplitude data

    Parameters
    ----------
    amp : array-like,
        amplitude of all 4 channels\n
    nbins : int, optional
        number of bins, 0 < nbins <= parameters.adcMax\n
    singlech : boolean, optional
        indicates whether the spectrum is single-channeled\n
    specRange : list, optional
        list for upper and lower bounds for the spectrum, in the form of `[lower, upper]` with 0 < lower < upper, or \
            `[[lower, upper], [lower, upper], [lower, upper], [lower, upper]]` spectrum ranges are channel-distinct\n
    specEdges : list or array-like or None, optional
        pre-defined spectrum edges, if not None, `specEdges` will be used instead of `nbins` and `specRange`
    binWidth : int or None, optional
        bin width to be used in the spectrum, within [1, parameters.adcMax], or None if not used\n

    Returns
    ----------
    spectrum : array-like or list of arrays,
        the corresponding spectrum(s) of the input\n
    x : array-like or list of arrays,
        bin centers of the spectrum(s)\n

    Raises
    ----------
    Exception
        when nbins is not within range `[1, parameters.adcMax]`, or specRange is not in correct form\n
    """

    sameRange = True
    if specEdges is None:
        if binWidth is None:
            if nbins <= 0 or nbins > parameters.adcMax:
                raise Exception('getSpectrum: parameter \'nbins\' out of range [1-' + str(parameters.adcMax) + ']')
        else:
            if binWidth <= 0 or binWidth > parameters.adcMax:
                raise Exception('getSpectrum: parameter \'binWidth\' out of range [1-' + str(parameters.adcMax) + ']')
        if isinstance(specRange[0], list):
            sameRange = False
            if singlech:
                raise Exception('getSpectrum: parameter \'specRange\' should not be channel-distinct for single-channeled data. Please specify only spectrum '
                    'range of the corresponding channel')
            if len(specRange) != 4:
                raise Exception('getSpectrum: when specifiying channel-distinct spectrum ranges, \'specRange\' should contain spectrum ranges for all 4 '
                    'channels')
            for ich in range(4):
                if len(specRange[ich]) != 2 or specRange[ich][0] < 0 or specRange[ich][1] <= specRange[ich][0]:
                    raise Exception('getSpectrum: parameter \'specRange\' for channel ' + str(ich) + ' is not given correct form')
        else:
            if len(specRange) != 2 or specRange[0] < 0 or specRange[1] <= specRange[0]:
                raise Exception('getSpectrum: parameter \'specRange\' is not given correct form')
    
    if not singlech:
        spectrum = []
        x = []
        for ich in range(4):
            if specEdges is None:
                curRange = (specRange[0], specRange[1]) if sameRange else (specRange[ich][0], specRange[ich][1])
                if binWidth is None:
                    specch, xch = np.histogram(amp[ich], bins = nbins, range = curRange)
                else:
                    binEdges = np.arange(curRange[0], curRange[1] + binWidth, binWidth)
                    specch, xch = np.histogram(amp[ich], bins = binEdges)
            else:
                specch, xch = np.histogram(amp[ich], bins = specEdges)
            spectrum.append(specch)
            x.append((xch[:-1] + xch[1:]) / 2)
    else:
        if specEdges is None:
            if binWidth is None:
                spectrum, x = np.histogram(amp, bins = nbins, range = (specRange[0], specRange[1]))
            else:
                binEdges = np.arange(specRange[0], specRange[1] + binWidth, binWidth)
                spectrum, x = np.histogram(amp, bins = binEdges)
        else:
            spectrum, x = np.histogram(amp, bins = specEdges)
        x = (x[:-1] + x[1:]) / 2
    return spectrum, x

def findPackPos(data, pattern, getLen = False):

    """
    Function for getting data pack positions according to given patterns

    Parameters
    ----------
    data : bytes or str,
        raw data\n
    pattern : Pattern,
        pattern of given data pack\n
    getLen : boolean, optional
        indicates whether the lengths of every pack will be returned

    Returns
    ----------
    pos : array-like,
        positions of data packs\n
    packLen : array-like, optional
        lengths of data packs
    """

    pos = []
    packLen = []
    for ip in pattern.finditer(data):
        pos.append(ip.start())
        packLen.append(ip.end() - ip.start())
    if getLen:
        return np.array(pos), np.array(packLen)
    else:
        return np.array(pos)

def extractFilename(filename):

    """
    Auxiliary function to extract file name from whole path, to make sure this programme works on both Windows and Unix-based systems

    Parameters
    ----------
    filename : str,
        file name with possible path suffix using either ' \ ' or ' / ' as seperator\n

    Returns
    ----------
    str,
        file name without path suffix
    """

    return re.split('\\\|/', filename)[-1]

def readECconfig(configFile, dataType = ''):

    """
    Function to read EC configurations from specified config file

    Parameters
    ----------
    configFile : str,
        name of the config file (with full path), namely `xrayConfig`, `sourceConfig` and `HPGeConfig` in gridParameters.py\n
    dataType : str,
        indicates type of data, with 'xray' being x-ray data, 'source' being source data and `hpge` being HPGe data

    Returns
    ----------
    fileList : list of dicts or dict,
        list of files to readout for all 4 channels for x-ray data or dict for source data and HPGe data, with list of each channel being a dictionary \
            in the form: \n
        files = {
                'date' :       list of experiment dates, for x-ray and HPGe data only
                'energy' :    list of corresponding energy, for x-ray and HPGe data only
                'source' :    list of corresponding sources, for source data only
                'data' :        list of data file names,
                'bkg' :         list of corresponding background data files, for x-ray and source data only
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
                configDict[channel] = {
                    'date' :             date, 
                    'energy' :          energy,
                    'data' :             dataFile,
                    'bkg' :               bkgFile,
                }

                
        #Source EC readout config
        elif dataType == 'source':
            source = []
            dataFile = []
            bkgFile = []
            for ifile in range(len(config)):
                filePath = config[ifile][0]
                source.append(config[ifile][1])
                dataFile.append(filePath + config[ifile][2])
                bkgFile.append(filePath + config[ifile][3])
            configDict = {
                'source' :          source,
                'data' :             dataFile,
                'bkg' :               bkgFile,
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
        name of the config file (with full path), namely `xrayProcessConfig` and `sourceProcessConfig` in gridParameters.py\n
    
    Returns
    ----------
    processOptions : tuple, 
        in the form of `(isHex, isCi, newProgramme, timeCut, rateStyle, deduplicate, doCorr, odr, maxiter, \
            bound)`\n
    the options include: \n
        isHex : boolean
            indicates whether the input file is in hexprint output format\n
        isCi : int
            indicates whether the input file has charge injection(CI) part, 0 for no CI, 1 for CI, 2 for multiple CI\n
        newProgramme : boolean
            indicates whether the data is produced with new hardware programme(6th ver.)\n
        timeCut : int
            cut of time data in seconds, specially designed for temp-bias data with pid bias control(6th ver.)\n
        rateStyle : str
            the style of count rate correction, '' for no corraction, 's' for correction with small data packs(512byte), 'l' for calculation \
                with large data packs(4096byte)\n
        deduplicate : boolean
            indicates whether the data will be deduplicated, for BINARY(.dat) data of new programme (6th ver.) only\n
        doCorr : boolean
            indicates whether the temperature-bias correction will be done, to avoid warning info in the output\n
        odr : boolean
            indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
        maxiter : int
            maximum number of iterations, 0 for auto range correction off\n
        bound : float
            boundary for auto rangecorrection in `\sigma`, with upper and lower bounds being `\mu - bound * \sigma` and `\mu + bound * \sigma`\n
    
    Raises
    ----------
    Exception
        when error reading config file (e.g. config file with wrong format)
    """

    isHex = False
    isCi = 0
    newProgramme = False
    timeCut = 0.
    rateStyle = ''
    deduplicate = False
    doCorr = False
    odr = False
    maxiter = 1
    bound = 3.

    try:
        with open(configFile, 'r') as fin:
            config = json.load(fin)
        #Determine whether the individual options are in correct type
        if isinstance(config['isHex'], bool):
            isHex = config['isHex']
        else:
            raise Exception()
        if isinstance(config['isCi'], int) and config['isCi'] < 3 and config['isCi'] > -1:
            isCi = config['isCi']
        else:
            raise Exception()
        if isinstance(config['newProgramme'], bool):
            newProgramme = config['newProgramme']
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
    except:
        raise Exception('readECProcessConfig: error reading config file \"' + configFile + '\"')

    return isHex, isCi, newProgramme, timeCut, rateStyle, deduplicate, doCorr, odr, maxiter, bound

def convertStrToTimestamp(timestr, format = '%Y%m%d', deltaTime = (0, 0, 0, 0)):

    """
    Auxiliary function to convert time string to timestamp according to specific format with specified time difference

    Parameters
    ----------
    timestr : str, 
        time string to convert\n
    format : str, optional
        corresponding format of the time string
    deltaTime : tuple. optional\n
        time difference to add to the UTC time
    
    Returns
    ----------
    timestamp : float, 
        timestamp corresponding to the input time string
    """

    currtime = datetime.datetime.strptime(timestr, format)
    currtime += datetime.timedelta(days = deltaTime[0], hours = deltaTime[1], minutes = deltaTime[2], seconds = deltaTime[3])
    return currtime.timestamp()

def convertTimestampToStr(timestamp, format = '%Y/%m/%d %H:%M:%S', deltaTime = (0, 0, 0, 0)):

    """
    Auxiliary function to convert timestamp to time string in specific format with specified time difference

    Parameters
    ----------
    timestamp : float, 
        timestamp to convert\n
    format : str, optional
        format of the time string to generate
    deltaTime : tuple. optional\n
        time difference to add to the UTC time
    
    Returns
    ----------
    timestr : str, 
        time string corresponding to the input timestamp
    """

    currtime = datetime.datetime.utcfromtimestamp(timestamp)
    currtime += datetime.timedelta(days = deltaTime[0], hours = deltaTime[1], minutes = deltaTime[2], seconds = deltaTime[3])
    return currtime.strftime(format)

def addDictToExistingData(newSciData, sciData, newTelData, telData, newScanData = {}, scanData = {}, isCi = 0, isScan = False):

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
    newScanData : dict, optional
        dictionary of new readout I-V scan data\n
    scanData : dict, optional
        dictionary of existing I-V scan data, possibly missing essential keys\n
    isCi : int, optional
        indicates whether the input file has charge injection(CI) part, 0 for no CI, 1 for CI, 2 for multiple CI\n
    isScan : boolean, optional
        indicates whether the input file has I-V scan part\n
    
    Returns
    ----------
    sciData : dict
        dictionary of existing science data with new data added\n
    telData : dict
        dictionary of existing telemetry data with new data added\n
    scanData : dict, optional
        dictionary of existing I-V scan data with new data added, or empty dict if `isScan` is `False`\n
    """

    #Process science data
    sciMultiChannel = ['amp', 'timestampEvt', 'sciNum', 'ampCI', 'timestampEvtCI']
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
    telMultiChannel = ['tempSipm', 'tempAdc', 'vMon', 'iMon', 'bias']
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
    
    #Process I-V scan data
    scanMultiChannel = ['vScan', 'iScan']
    newData = False
    if isScan:
        for key in newScanData:
            if key not in scanData:
                scanData.update({
                    key:        newScanData[key], 
                })
                if isCi == 2:
                    scanData.update({
                        'multiScan':        True, 
                    })
            else:
                newData = True
                if key in scanMultiChannel:
                    if not scanData['multiScan']:
                        newKeyData = []
                        for ich in range(4):
                            newKeyData.append([scanData[key][ich]])
                        scanData[key] = newKeyData
                    scanData[key] = list(scanData[key])
                    for ich in range(4):
                        scanData[key][ich] = list(scanData[key][ich])
                        if isCi == 2:
                            for isc in range(len(newScanData[key][ich])):
                                scanData[key][ich].append(newScanData[key][ich][isc])
                        else:
                            scanData[key][ich].append(newScanData[key][ich])
                    scanData[key] = np.array(scanData[key], dtype = object)
                else:
                    if not scanData['multiScan']:
                        scanData[key] = [scanData[key]]
                    scanData[key] = list(scanData[key])
                    if isCi == 2:
                        for isc in range(len(newScanData[key])):
                            scanData[key].append(newScanData[key][isc])
                    else:
                        scanData[key].append(newScanData[key])
                    scanData[key] = np.array(scanData[key], dtype = object)
        if newData:
            scanData.update({
                'multiScan':        True, 
            })
    
    #Output
    return sciData, telData, scanData

#*************************************************************************************************************************************************************
#**********************************************************************Basic I/O part**********************************************************************
#*************************************************************************************************************************************************************

def getTime(filename):

    """
    NOTE: CURRENTLY UNUSED
    ----------
    Function for reading out the total time from the filename\n

    Parameters
    ----------
    filename : str,
        name of the input file\n

    Returns
    ----------
    time : int,
        total time in seconds, or -1 if the filename cannot be parsed\n
    """

    time = 0
    lines = filename.split('_')
    for line in lines:
        if line.endswith('s'):
            #*** seconds
            try:
                time = int(line[:-1])
                return time
            except:
                pass
        elif line.endswith('m'):
            #*** minutes
            try:
                time = int(line[:-1]) * 60
                return time
            except:
                pass
        elif line.endswith('h'):
            #*** hours
            try:
                time = int(line[:-1]) * 3600
                return time
            except:
                pass
    print('getTime: unable to parse filename ' + filename)
    return -1

def crcCheck(data, crc):

    """
    Function for checking crc

    Parameters
    ----------
    data : list or array-like or bytes,
        data parts for calculating crc\n
    crc : list or array-like or bytes,
        crc parts in original data\n

    Returns
    ----------
    boolean,
        True if crc check is correct, False if not\n
    """

    lineBytes = b''
    if not isinstance(data, bytes):
        try:
            #Convert original data to bytes
            for x in data:
                lineBytes += struct.pack('B', x)
        except:
            return False
        crcCalculated = crc16.crc16xmodem(lineBytes)
    else:
        crcCalculated = crc16.crc16xmodem(data)
    crcData = crc[0] * 256 + crc[1]
    if not crcCalculated == crcData:
        return False
    return True

def deduplicate(sciPackData, telPackData = None, generateRaw = False):

    """
    NOTE: CURRENTLY FOR BINARY(.dat) FILES ONLY
    ----------
    Deduplication function for raw data (6th ver.)

    Parameters
    ----------
    sciPackData : list of bytes or list of str,
        raw science data packs to deduplicate\n
    telPackData : list of bytes or list of str or None, optional
        raw telemetry data packs to deduplicate, or None for single list deduplication\n
    generateRaw : boolean, optional
        indicates whether the deduplicated raw data will be recreated, if `False`, only the array of packs will be returned\n

    Returns
    ----------
    if generateRaw is `True`: \n
    deduplicated : bytes,
        deduplicated raw data\n
    else: \n
    sciPackData : array of bytes or list of str,
        deduplicated array of science data packs\n
    telPackData : array of bytes or list of str, optional
        deduplicated array of telemetry data packs\n
    """

    print('deduplicate: Deduplicating raw data...')

    #Deduplicate, with the original order of the data packs kept
    _, sciindex = np.unique(sciPackData, return_index = True)
    #Use dtype = object to avoid loss of '\x00' at the end of each bytes
    sciPackData = np.array(sciPackData, dtype = object)[np.sort(sciindex)]
    sciLen = len(sciPackData[0])
    if telPackData is not None:
        _, telindex = np.unique(telPackData, return_index = True)
        telPackData = np.array(telPackData, dtype = object)[np.sort(telindex)]
        telLen = len(telPackData[0])
    if generateRaw:
        #Generate deduplicated raw data
        print('Generating deduplicated raw data...')
        if telPackData is not None:
            deduplicated = bytearray(len(sciPackData) * sciLen + len(telPackData) * telLen)
        else:
            deduplicated = bytearray(len(sciPackData) * sciLen)
        #Telemetry data is generated first since telemetry packs are usually shorter in length
        currPos = 0
        if telPackData is not None:
            count = 0
            print('Telemetry data: ')
            for idata in telPackData:
                count += 1
                if parameters.deduplicateReportStep > 0 and count % parameters.deduplicateReportStep == 0:
                    #Status report
                    print(str(count) + '/' + str(len(telPackData)))
                deduplicated[currPos:currPos + telLen] = idata
                currPos += telLen
        print('Science data: ')
        count = 0
        for idata in sciPackData:
            count += 1
            if parameters.deduplicateReportStep > 0 and count % parameters.deduplicateReportStep == 0:
                #Status report
                print(str(count) + '/' + str(len(sciPackData)))
            deduplicated[currPos:currPos + sciLen] = idata
            currPos += sciLen
        return bytes(deduplicated)
    else:
        if telPackData is not None:
            return sciPackData, telPackData
        else:
            return sciPackData

def dataReadout(filename, isHex = False, isCi = 0, isScan = False, scanRange = [], rateStyle = '', newProgramme = False, timeCut = None, \
    isBin = True, doDeduplicate = False, reqNum = [], ignoreCrc = False, energyCut = None):

    """
    Function for reading out single Grid raw outout file, including: \n
    a) Single run/scan files, both text(.txt, both 5th. and 6th. ver, with hexprint supported) and binary(.dat 6th ver. only) files supported\n
    b) Multiple scans files(ascending timeline), currently only supporting text(.txt) files\n
    c) Multiple runs files(non-constant ascending timeline, or stacked single run files), currently only supporting binary(.dat) files\n

    Parameters
    ----------
    filename : str,
        name of the output file, currently supporting both text(.txt) and binary(.dat) files\n
    isHex : boolean, optional
        indicates whether the input file is in hexprint output format\n
    isCi : int, optional
        indicates whether the input file has charge injection(CI) part, 0 for no CI, 1 for CI, 2 for multiple CI\n
    isScan : boolean, optional
        indicates whether the input file has I-V scan part\n
    scanRange : list, optional
        containing the scan range for multiple scans, in the form of [firstScan, lastScan]\n
    rateStyle : str, optional
        the style of count rate correction, '' for no corraction, 's' for correction with small data packs(512byte), 'l' for calculation with large data \
            packs(4096byte)\n
    newProgramme : boolean, optional
        indicates whether the data is produced with new hardware programme(6th ver.)\n
    timeCut : int or list, optional
        cut of time data in seconds, specially designed for temp-bias data with pid bias control(6th ver.) or in-orbit data\n
    isBin: boolean, optional
        indicates whether the data is in binary form('*.dat' files)\n
    doDeduplicate : boolean, optional
        indicates whether the data will be deduplicated, for BINARY(.dat) data of new programme (6th ver.) only\n
    reqNum: list, optional
        required run numbers for science and telemetry data, in the form of [telreq, scireq]\n
    ignoreCrc : boolean, optional
        indicates whether CRC check will be ignored for current file. USE THIS ONLY FOR OLD (5TH VER.) NONE-HEXPRINT FILES\n
    energyCut : int or list, optional
        cut of energy in keV, specially designed for some in-orbit data\n

    Returns
    ----------
    sciExtracted : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciExtracted = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event uscount values, which corresponds to all events recorded in amp,
                    'timeCorrect':           array of arrays, correct live time, calculated with event uscount values, if rate correction is reqiured, 
                    'effectiveCount':       array of arrays, extracted effective count values, for new programme files only, 
                    'missingCount':          array of arrays, extracted missing count values, for new programme files only, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 

                    for data with charge injection (CI) part:
                    'ampCI':                    array of arrays, extracted amplitude data of CI part, 
                    'timestampEvtCI':      array of arrays, extracted event uscount values of CI part, which corresponds to all events recorded in ampCI, 
                    'effectiveCountCI':     array of arrays, extracted effective count values of CI part, for new programme files only, 
                    'missingCountCI':        array of arrays, extracted missing count values of CI part, for new programme files only, 
                }
        All data are in the form of `array(array, array, array, array)` for single scan(including multiple runs), or `array(array(array, len = NumOfScans), \
            array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans))` for multiple scans
    telExtracted : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telExtracted = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'tempAdc':                array of arrays, extracted ADC temperature, 
                    'vMon':                     array of arrays, extracted monitored voltage, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values(bias = vMon - iMon * R), 
                    'timestamp':             array-like or array of arrays, extracted uscount values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', \
                        for binary files only, 
                }
        Channel-seperated data are in the form of `array(array, array, array, array)` for single scan(including multiple runs), or \
            `array(array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans))` \
            for multiple scans\n
        Shared data are in the form of `array` for single scan(including multiple runs), or \
            `array(array, len = NumOfScans)` for multiple scans\n
    scanExtracted : dict, optional
        I-V scan data extracted from data packs, in the form of a dictionary as\n
        telExtracted = {
                    'vSet':                      array-like, set voltage, 
                    'vScan':                    array of arrays, monitored I-V scan voltages, in the form of array(array, array, array, array), 
                    'iScan':                     array of arrays, monitored I-V scan currents, in the form of array(array, array, array, array), 
                }

    Raises
    ----------
    Exception
        when given rateStyle is not in available styles\n
    """

    #------------------------------------------------------------Basic settings------------------------------------------------------------
    #Parse rate correction style
    styleAvailable = ['s', 'p', '']
    if not rateStyle in styleAvailable:
        raise Exception('dataReadout: count rate calculation style \'' + rateStyle + '\' not available')

    #Data corrrection for CRC check
    #Due to hardware problems, telemetry data packs in some files had values overwritten by CRC values, making standard
    #CRC checks impossible, only by searching for original values in previous and upcoming data packs can the CRC check be
    #done
    dataBuffer = [] #to fix the problem of telemetry data CRC check, checking the data packs before the current pack
    lineBuffer = [] #also for the DAMN telemetry data CRC check, checking the data packs after the current pack
    bufferLen = 500 #Search length for data packs before the current pack
    lineLen = 500 #Search length for data packs after the current pack

    print('dataReadout: processing \"' + extractFilename(filename) + '\"')

    #Read raw data from file
    if filename.endswith('.txt'):
        isBin = False
    if isBin:
        with open(filename, 'rb') as fin:
            lines = fin.read()
    else:
        with open(filename) as fin:
            lines = [line.rstrip() for line in fin]

    #Basic output
    amp = [] #after CI for data with CI
    tempSipm = []
    tempAdc = []
    vMon = []
    iMon = []
    bias = []
    utc = []
    uscount = [] #uscount in telemetry data
    uscountEvt = []
    timeCorrect = []
    effectiveCount = []
    missingCount = []
    sciNum = []
    telNum = []
    #CI output
    ampCI = []
    uscountEvtCI = []
    effectiveCountCI = []
    missingCountCI = []
    #IV scan output
    vSet = [] #i-v scan data
    vScan = [] #i-v scan data
    iScan = [] #i-v scan data

    #Set accepted scan range for multiple scans data
    ranged = False
    lower = -1
    upper = -1
    if isCi == 2 and len(scanRange) > 0:
        ranged = True
        lower = scanRange[0] - 1
        upper = scanRange[1] - 1

    #Science scan data, in all 4 channels
    for ich in range(4):
        amp.append([])
        ampCI.append([])
        tempSipm.append([])
        tempAdc.append([])
        vMon.append([])
        iMon.append([])
        bias.append([])
        uscountEvt.append([])
        uscountEvtCI.append([])
        vScan.append([])
        iScan.append([])
        sciNum.append([])

    #Basic settings
    bCi = True
    if isCi == 0:
        #Single scan
        bCi = False
    if isBin:
        #Binary files
        isHex = True
    if isHex:
        #Hexprint files, which have no CI and I-V scan part
        isCi = 0
        isScan = False
    nScan = -1
    indexOut = 0 #count of events with channel index out of range[1-4]
    crcError = 0 #count of crc error data
    utcError = 0 #count of utc error data
    lastuscount = -10 #latest telemetry uscount value in ascending sequence, used to split telemetry runs
    scilastuscount = -10 #latest event uscount value in ascending sequence, used to split science runs
    section = 0 #telemetry run number(section)
    scisection = 0 #science run number(section)
    reqRead = False #if required run specified, readout ends when required run is read
    #Beginning and ends(in uscount values) of science and telemetry runs
    sciBegin = []
    sciEnd = []
    telBegin = []
    telEnd = []

    #Set required run for binary files
    telreq = -1
    scireq = -1
    if isBin:
        if len(reqNum) == 2:
            telreq = reqNum[0]
            scireq = reqNum[1]
        elif extractFilename(filename) in parameters.reqFileRef:
            telreq = parameters.reqFileRef[extractFilename(filename)][0]
            scireq = parameters.reqFileRef[extractFilename(filename)][1]
        else:
            telreq = -1
            scireq = -1


    #------------------------------------------------------------Readout for binary files------------------------------------------------------------
    if isBin:
        #Patterns of science and telemetry data packs
        sciPattern = re.compile(parameters.pattrenRef['binSci'], re.S)
        telPattern = re.compile(parameters.pattrenRef['binTel'], re.S)
        #Find positions of all science and telemetry data packs
        sciPackPos = findPackPos(lines, sciPattern)
        telPackPos = findPackPos(lines, telPattern)
        #Get all science and telemetry data packs
        sciPackData = []
        telPackData = []
        for idata in range(len(sciPackPos)):
            sciPackData.append(lines[sciPackPos[idata]:sciPackPos[idata] + 512])
        for idata in range(len(telPackPos)):
            telPackData.append(lines[telPackPos[idata]:telPackPos[idata] + 512])
        #Do deduplicate data and generate deduplicated raw data, without readout
        if parameters.deduplicateFile:
            lines = deduplicate(sciPackData, telPackData, parameters.deduplicateFile)
            filenameNoPath = extractFilename(filename)
            filenameNoAppend = filenameNoPath[0:filenameNoPath.find(filenameNoPath.split('.')[-1]) - 1]
            with open(filenameNoAppend + '_deduplicated.' + filenameNoPath.split('.')[-1], 'wb+') as fout:
                fout.write(lines)
        else:
            if doDeduplicate:
                #Deduplicate data
                sciPackData, telPackData = deduplicate(sciPackData, telPackData)

            utcLaunch = 0.
            if parameters.utcCheck:
                utcLaunch = convertStrToTimestamp(parameters.launchTime)
            telcuruscount = 0
            #Readout of telemetry data packs
            for telData in telPackData:
                if len(telData) < 498:
                    #Avoid error caused by incomplete data packs
                    crcError += 1#49140 wrt this line, 53173 with this line, actual 49341
                    continue
                #CRC check
                if not crcCheck(telData[:496], telData[496:498]):
                    crcError += 1
                    continue
                for it in range(7):
                    #Get current uscount value
                    telcuruscount = struct.unpack('>Q', telData[15 + 70 * it:23 + 70 * it])[0] / parameters.internalFreq
                    if lastuscount < 0:
                        #Add first telemetry run begin time
                        telBegin.append(telcuruscount)
                    if telcuruscount < lastuscount:
                        #Non-ascending uscount value, which indicates run split point
                        section += 1
                        if telreq > -1 and section > telreq and scireq > -1 and scisection > scireq:
                            #Exit when the required run is read
                            reqRead = True
                        print('Telemetry scan #' + str(section))
                        telBegin.append(telcuruscount)
                        if lastuscount > -0.5:
                            #Add end time of current telemetry run
                            telEnd.append(lastuscount)
                    lastuscount = telcuruscount
                    curutc = struct.unpack('>L', telData[3 + 70 * it:7 + 70 * it])[0]
                    if parameters.utcCheck:
                        #UTC check for in-orbit data
                        if curutc <= utcLaunch:
                            utcError += 1
                            continue
                    if telreq < 0 or section == telreq:
                        #Extract data, only when no telemetry run is required or current run is the required run
                        telNum.append(section)
                        utc.append(curutc)
                        uscount.append(telcuruscount)
                        for ich in range(4):
                            curtemp = float(struct.unpack('>H', telData[23 + 2 * ich + 70 * it:25 + 2 * ich + 70 * it])[0])
                            tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                            curtemp = float(struct.unpack('>H', telData[31 + 2 * ich + 70 * it:33 + 2 * ich + 70 * it])[0])
                            tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                            vMon[ich].append(float(struct.unpack('>H', telData[39 + 2 * ich + 70 * it:41 + 2 * ich + 70 * it])[0]) / 4096.0 * 3.3 * 11.0)
                            iMon[ich].append(float(struct.unpack('>H', telData[47 + 2 * ich + 70 * it:49 + 2 * ich + 70 * it])[0]) / 4096.0 * 3.3)
                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * parameters.internalResistance)
                #Exit when the required run is read
                if reqRead:
                    break

            scicuruscount = 0
            #Readout of science data packs
            for sciData in sciPackData:
                if len(sciData) < 512:
                    #Avoid error caused by incomplete data packs
                    crcError += 1
                    continue
                timeEvtBegin = 0.0
                timeEvtEnd = 0.0
                timeIntv = []
                #CRC check
                if not crcCheck(sciData[:510], sciData[510:512]):
                    crcError += 1
                    continue
                #Channel of currrent event
                ch = sciData[3]
                if not newProgramme:
                    ch -= 1
                #Get first uscount of current package
                scicuruscount = float(struct.unpack('>Q', sciData[4:12])[0]) / parameters.internalFreq
                if ch > -1 and ch < 4:
                    if scilastuscount < 0:
                        #Add first science run begin time
                        sciBegin.append(scicuruscount)
                    if scicuruscount < scilastuscount:
                        #Non-ascending uscount value, which indicates run split point
                        scisection += 1
                        if telreq > -1 and section > telreq and scireq > -1 and scisection > scireq:
                            #Exit when the required run is read
                            reqRead = True
                        print('Science scan #' + str(scisection))
                        sciBegin.append(scicuruscount)
                        if scilastuscount > -0.5:
                            #Add end time of current science run
                            sciEnd.append(scilastuscount)
                    scilastuscount = scicuruscount
                    if scireq < 0 or scisection == scireq:
                        #Extract data, only when no sicence run is required or current run is the required run
                        sciNum[ch].append(scisection)
                        amp[ch].append(struct.unpack('>H', sciData[12:14])[0])
                        uscountEvt[ch].append(scicuruscount)
                else:
                    #Channel index out of bound [0-3]
                    indexOut += 1
                if rateStyle == 's':
                    #512B package rate correction, record beginning of all 512B data packs
                    timeIntv.append(scicuruscount)
                elif rateStyle == 'p':
                    #4096B package rate correction, record beginning of all 4096B data packs
                    timeEvtBegin = scicuruscount
                for ie in range(43):
                    #Channel of current event
                    ch = sciData[26 + 11 * ie]
                    if not newProgramme:
                        ch -= 1
                    #Get uscount of current event
                    scicuruscount = float(struct.unpack('>Q', sciData[27 + 11 * ie:35 + 11 * ie])[0]) / parameters.internalFreq
                    if ch > -1 and ch < 4:
                        if scicuruscount < scilastuscount:
                            #Non-ascending uscount value, which indicates run split point
                            scisection += 1
                            if telreq > -1 and section > telreq and scireq > -1 and scisection > scireq:
                                #Exit when the required run is read
                                reqRead = True
                            print('Science scan #' + str(scisection))
                            sciBegin.append(scicuruscount)
                            if scilastuscount > -0.5:
                                #Add end time of current science run
                                sciEnd.append(scilastuscount)
                        scilastuscount = scicuruscount
                        if scireq < 0 or scisection == scireq:
                            #Extract data, only when no sicence run is required or current run is the required run
                            sciNum[ch].append(scisection)
                            amp[ch].append(struct.unpack('>H', sciData[35 + 11 * ie:37 + 11 * ie])[0])
                            uscountEvt[ch].append(scicuruscount)
                    else:
                        #Channel index out of bound [0-3]
                        indexOut += 1
                    if rateStyle == 's':
                        #512B package rate correction, record beginning of all 512B data packs
                        timeIntv.append(scicuruscount)
                if rateStyle == 's':
                    #512B package rate correction, calculate correct live time of each 512B package
                    timeIntv = np.array(timeIntv)
                    timeCorrect += list(timeIntv[1:] - timeIntv[:-1])
                elif rateStyle == 'p':
                    #4096B package rate correction, calculate correct live time of each 4096B package
                    timeEvtEnd = scicuruscount
                    timeCorrect.append(timeEvtEnd - timeEvtBegin)
                if newProgramme:
                    effectiveCount.append(struct.unpack('>L', sciData[499:503])[0])
                    missingCount.append(struct.unpack('>L', sciData[503:507])[0])
                #Exit when the required run is read
                if reqRead:
                    break

    #------------------------------------------------------------Readout for text files------------------------------------------------------------
    else:
        #--------------------------------------------------Hexprint format files--------------------------------------------------
        if isHex:
            #Patterns of science and telemetry data packs
            if newProgramme:
                sciPattern = re.compile(parameters.pattrenRef['hexSciNew'], re.I)
                telPattern = re.compile(parameters.pattrenRef['hexTelNew'], re.I)
            else:
                sciPattern = re.compile(parameters.pattrenRef['hexSci'], re.I)
                telPattern = re.compile(parameters.pattrenRef['hexTel'], re.I)
            #Find positions of all science and telemetry data packs
            data = lines[0]
            sciPackPos = findPackPos(data, sciPattern)
            telPackPos = findPackPos(data, telPattern)
            #Get all science and telemetry data packs
            sciPackData = []
            telPackData = []
            for idata in range(len(sciPackPos)):
                sciPackData.append(data[sciPackPos[idata]:sciPackPos[idata] + 512 * 3])
            for idata in range(len(telPackPos)):
                telPackData.append(data[telPackPos[idata]:telPackPos[idata] + 512 * 3])
            if doDeduplicate:
                #Deduplicate data
                sciPackData, telPackData = deduplicate(sciPackData, telPackData)

            scicuruscount = 0
            #Readout of science data packs
            for sciData in sciPackData:
                #Extract integer sequence
                lineList = sciData.split(' ')
                lineFloat = []
                try:
                    for linestr in lineList:
                        lineFloat.append(int(linestr, 16))
                except:
                    pass
                if len(lineFloat) < 512:
                    continue
                timeEvtBegin = 0.0
                timeEvtEnd = 0.0
                timeIntv = []
                #CRC check for old programme(5th ver.)
                if not newProgramme:
                    #Remove first buffer when buffer reaches maximum size
                    if len(dataBuffer) >= bufferLen:
                        del dataBuffer[0]
                    #Fill buffer for checking upcoming data packs
                    dataBuffer.append(lineFloat[496:504])
                    #Check previous telemetry data
                    for buf in lineBuffer:
                        if buf[498:504] == lineFloat[498:504]:
                            #If a match is found, calculate correct CRC for the matched previous telemetry pack
                            bufcrc = buf[496:498]
                            buf[496:498] = lineFloat[496:498]
                            #Check CRC of the matched previous telemetry pack
                            if not crcCheck(buf[0:510], bufcrc):
                                #Remove the buffer if CRC check not passed
                                crcError += 1
                                del lineBuffer[lineBuffer.index(buf)]
                            else:
                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                for it in range(7):
                                    utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                    uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                    for ich in range(4):
                                        curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                        tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                        tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                        iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                        bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * parameters.internalResistance)
                                del lineBuffer[lineBuffer.index(buf)]
                            break
                    #CRC check of current science pack
                    if not crcCheck(lineFloat[:502], lineFloat[502:504]):
                        crcError += 1
                        continue
                else:
                    #CRC check of current science pack
                    if not crcCheck(lineFloat[:510], lineFloat[510:512]):
                        crcError += 1
                        continue
                #Channel of currrent event
                ch = lineFloat[3]
                #Get first uscount of current package
                scicuruscount = float(sum([lineFloat[4 + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq
                if newProgramme:
                    ch += 1
                if ch > 0 and ch < 5:
                    amp[ch - 1].append(lineFloat[12] * 256 + lineFloat[13])
                    uscountEvt[ch - 1].append(scicuruscount)
                else:
                    #Channel index out of bound [0-3]
                    indexOut += 1
                if rateStyle == 's':
                    #512B package rate correction, record beginning of all 512B data packs
                    timeIntv.append(scicuruscount)
                elif rateStyle == 'p':
                    #4096B package rate correction, calculate correct live time of each 4096B package
                    timeEvtBegin = scicuruscount
                for ie in range(43):
                    #Channel of currrent event
                    ch = lineFloat[26 + 11 * ie]
                    #Get uscount of current event
                    scicuruscount = float(sum([lineFloat[27 + 11 * ie + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq
                    if newProgramme:
                        ch += 1
                    if ch > 0 and ch < 5:
                        amp[ch - 1].append(lineFloat[35 + 11 * ie] * 256 + lineFloat[36 + 11 * ie])
                        uscountEvt[ch - 1].append(scicuruscount)
                    else:
                        indexOut += 1
                    if rateStyle == 's':
                        #512B package rate correction, record beginning of all 512B data packs
                        timeIntv.append(scicuruscount)
                if rateStyle == 's':
                    #512B package rate correction, calculate correct live time of each 512B package
                    timeIntv = np.array(timeIntv)
                    timeCorrect += list(timeIntv[1:] - timeIntv[:-1])
                elif rateStyle == 'p':
                    #4096B package rate correction, calculate correct live time of each 4096B package
                    timeEvtEnd = scicuruscount
                    timeCorrect.append(timeEvtEnd - timeEvtBegin)
                if newProgramme:
                    effectiveCount.append(int(sum([lineFloat[499 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                    missingCount.append(int(sum([lineFloat[503 + ic] * 256 ** (3 - ic) for ic in range(4)])))

            telcuruscount = 0
            #Readout of telemetry data packs
            for telData in telPackData:
                #Extract integer sequence
                lineList = telData.split(' ')
                lineFloat = []
                try:
                    for linestr in lineList:
                        lineFloat.append(int(linestr, 16))
                except:
                    pass
                if len(lineFloat) < 512:
                    continue
                #CRC check for old programme(5th ver.)
                if not newProgramme:
                    #Check the data before the current data for match to correct current data
                    crcCorrect = []
                    for i in dataBuffer:
                        if i[2:] == lineFloat[498:504]:
                            crcCorrect = i[:2]
                            break
                    #Check the previous telemetry data
                    for buf in lineBuffer:
                        if buf[498:504] == lineFloat[498:504]:
                            #If a match is found, calculate correct CRC for the matched previous telemetry pack
                            bufcrc = buf[496:498]
                            buf[496:498] = lineFloat[496:498]
                            #Check CRC of the matched previous telemetry pack
                            if not crcCheck(buf[0:510], bufcrc):
                                #Remove the buffer if CRC check not passed
                                crcError += 1
                                del lineBuffer[lineBuffer.index(buf)]
                            else:
                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                for it in range(7):
                                    utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                    uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                    for ich in range(4):
                                        curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                        tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                        tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                        iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                        bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * parameters.internalResistance)
                                del lineBuffer[lineBuffer.index(buf)]
                            break
                    #CRC check of current telemetry pack and buffer fill
                    #Remove first buffer when buffer reaches maximum size
                    if len(dataBuffer) >= bufferLen:
                        del dataBuffer[0]
                    #Fill buffer for checking upcoming data packs
                    dataBuffer.append(lineFloat[496:504])
                    #Store the CRC values of current telemetry pack
                    temp = lineFloat[496:498]
                    if not len(crcCorrect) == 0:
                        #If a match is found in the previous data pack, use it to correct current telemetry data pack
                        lineFloat[496:498] = crcCorrect
                    else:
                        #If no match is found, fill current telemetry pack to to-be-checked buffer
                        if len(lineBuffer) >= lineLen:
                            del lineBuffer[0]
                            crcError += 1
                        lineBuffer.append(lineFloat[:512])
                    #CRC check of current telemetry pack
                    if not crcCheck(lineFloat[:510], temp):
                        crcError += 1
                        continue
                else:
                    #CRC check of current telemetry pack
                    if not crcCheck(lineFloat[:496], lineFloat[496:498]):
                        crcError += 1
                        continue
                for it in range(7):
                    utc.append(float(sum([lineFloat[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                    uscount.append(float(sum([lineFloat[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                    for ich in range(4):
                        curtemp = float(lineFloat[23 + 2 * ich + 70 * it] * 256 + lineFloat[24 + 2 * ich + 70 * it])
                        tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                        curtemp = float(lineFloat[31 + 2 * ich + 70 * it] * 256 + lineFloat[32 + 2 * ich + 70 * it])
                        tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                        vMon[ich].append(float(lineFloat[39 + 2 * ich + 70 * it] * 256 + lineFloat[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                        iMon[ich].append(float(lineFloat[47 + 2 * ich + 70 * it] * 256 + lineFloat[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3 / 2.0)
                        bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * parameters.internalResistance)

        #--------------------------------------------------Non-hexprint(normal) format files--------------------------------------------------
        else:
            #For some files, ignore CRC check to ensure there is still telemetry data for further processing
            if extractFilename(filename) in parameters.ignoreCrcRef:
                ignoreCrc = True
            
            for line in lines:
                #Readout of I-V scan data
                if isScan and 'Point' in line:
                    lineList = line.split(',')
                    scanData = []
                    bScan = True #ensuring that corrupted data are disposed
                    try:
                        scanData.append(int(lineList[3]))
                        scanList = lineList[4].split()
                        for ich in range(4):
                            scanData.append(float(scanList[2 * ich]) / 4096.0 * 3.3 * 11.0)
                            scanData.append(float(scanList[1 + 2 * ich]) / 4096.0 * 2.0 * 3.3)
                    except:
                        bScan = False
                    if bScan:
                        vSet.append(scanData[0])
                        for ich in range(4):
                            vScan[ich].append(scanData[1 + 2 * ich])
                            iScan[ich].append([scanData[2 + 2 * ich]])
                    continue

                #Begin and end markers of CI, also scan counts for multiple scans
                if isCi == 2 and 'Begin' in line:
                    #For multiple scans, beginning markers of CI indicate beginning of a new scan
                    bCi = True
                    nScan += 1
                    if ranged and not (nScan >= lower and nScan <= upper):
                        #Skip scans outside the required scan range
                        print('Skipping scan #' + str(nScan + 1))
                    else:
                        print('Scan #' + str(nScan + 1))
                    uscount.append([])
                    utc.append([])
                    if not rateStyle == '':
                        timeCorrect.append([])
                    if newProgramme:
                        effectiveCount.append([])
                        effectiveCountCI.append([])
                        missingCount.append([])
                        missingCountCI.append([])
                    for ich in range(4):
                        amp[ich].append([])
                        ampCI[ich].append([])
                        tempSipm[ich].append([])
                        tempAdc[ich].append([])
                        vMon[ich].append([])
                        iMon[ich].append([])
                        bias[ich].append([])
                        uscountEvt[ich].append([])
                        uscountEvtCI[ich].append([])
                    continue
                elif 'End' in line:
                    bCi = False
                    continue

                #Readout of single line
                lineList = line.split(' ')
                if len(lineList) > 502:
                    #Non-hexprint(normal) format files
                    if isCi == 2 and ranged:
                        if not (nScan >= lower and nScan <= upper):
                            continue
                    lineFloat = []
                    try:
                        for linestr in lineList:
                            lineFloat.append(int(linestr))
                    except:
                        pass

                    for il in range(len(lineFloat)):
                        #Readout of science data pack
                        if (lineFloat[il] == 170 and lineFloat[il + 1] == 187 and lineFloat[il + 2] == 204) and \
                            ((il + 502 <= len(lineFloat) and lineFloat[il + 499] == 221 and lineFloat[il + 500] == 238 and lineFloat[il + 501] == 255 and \
                            (not newProgramme)) or (il + 510 <= len(lineFloat) and lineFloat[il + 507] == 221 and lineFloat[il + 508] == 238 and \
                            lineFloat[il + 509] == 255 and newProgramme)):
                                timeEvtBegin = 0.0
                                timeEvtEnd = 0.0
                                timeIntv = []
                                #CRC check for old programme(5th ver.)
                                if not newProgramme and not ignoreCrc:
                                    #CRC check and CRC buffer fill
                                    if len(dataBuffer) >= bufferLen:
                                        del dataBuffer[0]
                                    #Fill buffer for checking upcoming data packs
                                    dataBuffer.append(lineFloat[496:504])
                                    #Check previous telemetry data
                                    for buf in lineBuffer:
                                        if buf[498:504] == lineFloat[498:504]:
                                            #If a match is found, calculate correct CRC for the matched previous telemetry pack
                                            bufMatch = buf
                                            bufcrc = buf[496:498]
                                            bufMatch[496:498] = lineFloat[496:498]
                                            #Check CRC of the matched previous telemetry pack
                                            if not crcCheck(bufMatch[0:510], bufcrc):
                                                #Remove the buffer if CRC check not passed
                                                crcError += 1
                                                del lineBuffer[lineBuffer.index(buf)]
                                            else:
                                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                                if isCi == 2:
                                                    for it in range(7):
                                                        utc[nScan].append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount[nScan].append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                                        for ich in range(4):
                                                            curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                                            tempSipm[ich][nScan].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                                            tempAdc[ich][nScan].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            vMon[ich][nScan].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich][nScan].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich][nScan].append(vMon[ich][nScan][-1] - iMon[ich][nScan][-1] * parameters.internalResistance)
                                                else:
                                                    for it in range(7):
                                                        utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                                        for ich in range(4):
                                                            curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                                            tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                                            tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * parameters.internalResistance)
                                                del lineBuffer[lineBuffer.index(buf)]
                                            break
                                    #CRC check of current science pack
                                    if not crcCheck(lineFloat[:502], lineFloat[502:504]):
                                        crcError += 1
                                        continue
                                else:
                                    #CRC check of current science pack
                                    if not crcCheck(lineFloat[:510], lineFloat[510:512]):
                                        crcError += 1
                                        continue
                                #Channel of currrent event
                                ch = lineFloat[il + 3]
                                #Get first uscount of current package
                                scicuruscount = float(sum([lineFloat[il + 4 + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq
                                if newProgramme:
                                    ch += 1
                                if isCi == 2:
                                    if ch > 0 and ch < 5:
                                        if bCi:
                                            ampCI[ch - 1][nScan].append(lineFloat[il + 12] * 256 + lineFloat[il + 13])
                                            uscountEvtCI[ch - 1][nScan].append(scicuruscount)
                                        else:
                                            amp[ch - 1][nScan].append(lineFloat[il + 12] * 256 + lineFloat[il + 13])
                                            uscountEvt[ch - 1][nScan].append(scicuruscount)
                                    else:
                                        #Channel index out of bound [0-3]
                                        indexOut += 1
                                else:
                                    if ch > 0 and ch < 5:
                                        if bCi:
                                            ampCI[ch - 1].append(lineFloat[il + 12] * 256 + lineFloat[il + 13])
                                            uscountEvtCI[ch - 1].append(scicuruscount)
                                        else:
                                            amp[ch - 1].append(lineFloat[il + 12] * 256 + lineFloat[il + 13])
                                            uscountEvt[ch - 1].append(scicuruscount)
                                    else:
                                        #Channel index out of bound [0-3]
                                        indexOut += 1
                                if rateStyle == 's' and not bCi:
                                    #512B package rate correction, record beginning of all 512B data packs
                                    timeIntv.append(scicuruscount)
                                elif rateStyle == 'p' and not bCi:
                                    #4096B package rate correction, calculate correct live time of each 4096B package
                                    timeEvtBegin = scicuruscount
                                for ie in range(43):
                                    #Channel of currrent event
                                    ch = lineFloat[il + 26 + 11 * ie]
                                    #Get uscount of current event
                                    scicuruscount = float(sum([lineFloat[il + 27 + 11 * ie + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq
                                    if newProgramme:
                                        ch += 1
                                    if isCi == 2:
                                        if ch > 0 and ch < 5:
                                            if bCi:
                                                ampCI[ch - 1][nScan].append(lineFloat[il + 35 + 11 * ie] * 256 + lineFloat[il + 36 + 11 * ie])
                                                uscountEvtCI[ch - 1][nScan].append(scicuruscount)
                                            else:
                                                amp[ch - 1][nScan].append(lineFloat[il + 35 + 11 * ie] * 256 + lineFloat[il + 36 + 11 * ie])
                                                uscountEvt[ch - 1][nScan].append(scicuruscount)
                                        else:
                                            indexOut += 1
                                    else:
                                        if ch > 0 and ch < 5:
                                            if bCi:
                                                ampCI[ch - 1].append(lineFloat[il + 35 + 11 * ie] * 256 + lineFloat[il + 36 + 11 * ie])
                                                uscountEvtCI[ch - 1].append(scicuruscount)
                                            else:
                                                amp[ch - 1].append(lineFloat[il + 35 + 11 * ie] * 256 + lineFloat[il + 36 + 11 * ie])
                                                uscountEvt[ch - 1].append(scicuruscount)
                                        else:
                                            indexOut += 1
                                    if rateStyle == 's' and not bCi:
                                        #512B package rate correction, record beginning of all 512B data packs
                                        timeIntv.append(scicuruscount)
                                if rateStyle == 's' and not bCi:
                                    #512B package rate correction, calculate correct live time of each 512B package
                                    timeIntv = np.array(timeIntv)
                                    if isCi == 2:
                                        timeCorrect[nScan] += list(timeIntv[1:] - timeIntv[:-1])
                                    else:
                                        timeCorrect += list(timeIntv[1:] - timeIntv[:-1])
                                elif rateStyle == 'p' and not bCi:
                                    #4096B package rate correction, calculate correct live time of each 4096B package
                                    timeEvtEnd = scicuruscount
                                    if isCi == 2:
                                        timeCorrect[nScan].append(timeEvtEnd - timeEvtBegin)
                                    else:
                                        timeCorrect.append(timeEvtEnd - timeEvtBegin)
                                if newProgramme:
                                    if isCi == 2:
                                        if bCi:
                                            effectiveCountCI[nScan].append(int(sum([lineFloat[il + 499 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                                            missingCountCI[nScan].append(int(sum([lineFloat[il + 503 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                                        else:
                                            effectiveCount[nScan].append(int(sum([lineFloat[il + 499 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                                            missingCount[nScan].append(int(sum([lineFloat[il + 503 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                                    else:
                                        if bCi:
                                            effectiveCountCI.append(int(sum([lineFloat[il + 499 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                                            missingCountCI.append(int(sum([lineFloat[il + 503 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                                        else:
                                            effectiveCount.append(int(sum([lineFloat[il + 499 + ic] * 256 ** (3 - ic) for ic in range(4)])))
                                            missingCount.append(int(sum([lineFloat[il + 503 + ic] * 256 ** (3 - ic) for ic in range(4)])))

                        #Readout of telemetry data pack
                        elif (lineFloat[il] == 1 and lineFloat[il + 1] == 35 and lineFloat[il + 2] == 69 and il + 502 <= len(lineFloat) and \
                            lineFloat[il + 493] == 103 and lineFloat[il + 494] == 137 and lineFloat[il + 495] == 16 and (not newProgramme)) or \
                            (lineFloat[il] == 18 and lineFloat[il + 1] == 52 and lineFloat[il + 2] == 86 and il + 502 <= len(lineFloat) and \
                            lineFloat[il + 493] == 120 and lineFloat[il + 494] == 154 and lineFloat[il + 495] == 188 and newProgramme):
                                #CRC check for old programme(5th ver.)
                                if not newProgramme and not ignoreCrc:
                                    #Check the data before the current data for match to correct current data
                                    crcCorrect = []
                                    for databuf in dataBuffer:
                                        if databuf[2:] == lineFloat[498:504]:
                                            crcCorrect = databuf[:2]
                                            break
                                    #Check the previous telemetry data
                                    for buf in lineBuffer:
                                        if buf[498:504] == lineFloat[498:504]:
                                            #If a match is found, calculate correct CRC for the matched previous telemetry pack
                                            bufMatch = buf
                                            bufcrc = buf[496:498]
                                            bufMatch[496:498] = lineFloat[496:498]
                                            #Check CRC of the matched previous telemetry pack
                                            if not crcCheck(bufMatch[0:510], bufcrc):
                                                #Remove the buffer if CRC check not passed
                                                crcError += 1
                                                del lineBuffer[lineBuffer.index(buf)]
                                            else:
                                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                                if isCi == 2:
                                                    for it in range(7):
                                                        utc[nScan].append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount[nScan].append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                                        for ich in range(4):
                                                            tempSipm[ich][nScan].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] > 2048 \
                                                                else tempSipm[ich][nScan].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it]) / 16.0)
                                                            tempAdc[ich][nScan].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] > 2048 \
                                                                else tempAdc[ich][nScan].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it]) / 16.0)
                                                            vMon[ich][nScan].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich][nScan].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich][nScan].append(vMon[ich][nScan][-1] - iMon[ich][nScan][-1] * parameters.internalResistance)
                                                else:
                                                    for it in range(7):
                                                        utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                                        for ich in range(4):
                                                            tempSipm[ich].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] > 2048 \
                                                                else tempSipm[ich].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it]) / 16.0)
                                                            tempAdc[ich].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] > 2048 \
                                                                else tempAdc[ich].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it]) / 16.0)
                                                            vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * parameters.internalResistance)
                                                del lineBuffer[lineBuffer.index(buf)]
                                            break
                                    #CRC check of current telemetry pack and buffer fill
                                    #Remove first buffer when buffer reaches maximum size
                                    if len(dataBuffer) >= bufferLen:
                                        del dataBuffer[0]
                                    #Fill buffer for checking upcoming data packs
                                    dataBuffer.append(lineFloat[496:504])
                                    #Store the CRC values of current telemetry pack
                                    temp = lineFloat[496:498]
                                    if not len(crcCorrect) == 0:
                                        #If a match is found in the previous data pack, use it to correct current telemetry data pack
                                        lineFloat[496:498] = crcCorrect
                                    else:
                                        #If no match is found, fill current telemetry pack to to-be-checked buffer
                                        if len(lineBuffer) >= lineLen:
                                            del lineBuffer[0]
                                            crcError += 1
                                        lineBuffer.append(lineFloat)
                                    #CRC check of current telemetry pack
                                    if not crcCheck(lineFloat[0:510], temp):
                                        crcError += 1
                                        continue
                                else:
                                    #CRC check of current telemetry pack
                                    if not crcCheck(lineFloat[0:496], lineFloat[496:498]) and not ignoreCrc:
                                        crcError += 1
                                        continue
                                for it in range(7):
                                    if isCi == 2:
                                        utc[nScan].append(float(sum([lineFloat[il + 3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                        uscount[nScan].append(float(sum([lineFloat[il + 15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                        for ich in range(4):
                                            tempSipm[ich][nScan].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] > 2048 \
                                                else tempSipm[ich][nScan].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it]) / 16.0)
                                            tempAdc[ich][nScan].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] > 2048 \
                                                else tempAdc[ich][nScan].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it]) / 16.0)
                                            vMon[ich][nScan].append(float(lineFloat[il + 39 + 2 * ich + 70 * it] * 256 + lineFloat[il + 40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                            iMon[ich][nScan].append(float(lineFloat[il + 47 + 2 * ich + 70 * it] * 256 + lineFloat[il + 48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                            bias[ich][nScan].append(vMon[ich][nScan][-1] - iMon[ich][nScan][-1] * parameters.internalResistance)
                                    else:
                                        utc.append(float(sum([lineFloat[il + 3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                        uscount.append(float(sum([lineFloat[il + 15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / parameters.internalFreq)
                                        for ich in range(4):
                                            tempSipm[ich].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] > 2048 \
                                                else tempSipm[ich].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it]) / 16.0)
                                            tempAdc[ich].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] > 2048 \
                                                else tempAdc[ich].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it]) / 16.0)
                                            vMon[ich].append(float(lineFloat[il + 39 + 2 * ich + 70 * it] * 256 + lineFloat[il + 40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                            iMon[ich].append(float(lineFloat[il + 47 + 2 * ich + 70 * it] * 256 + lineFloat[il + 48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * parameters.internalResistance)


    #------------------------------------------------------------Readout data preliminary processing------------------------------------------------------------
    if not parameters.deduplicateFile:
        errorCount = crcError + len(lineBuffer)
        print(str(errorCount) + ' data packs with crc error')
        print(str(indexOut) + ' events with channel out of bound[0-3]')
        if parameters.utcCheck:
            print(str(utcError) + ' data packs with utc error')
        
        if isBin:
            if len(uscount) == 0:
                print('No telemetry data was read. ')
                if telreq > -1 and section <= telreq:
                    print('    Due to telemetry run number ' + str(telreq) + ' exceeding total run number of ' + str(section))
            for ich in range(4):
                if len(uscountEvt[ich]) == 0:
                    print('No science data was read for channel ' + str(ich) + '. ')
                    if scireq > -1 and scisection <= scireq:
                        print('    Due to science run number ' + str(scireq) + ' exceeding total run number of ' + str(scisection))
            telEnd.append(uscount[-1])
            sciEnd.append(np.max([np.max(uscountEvt[ich]) for ich in range(4)]))

            #Display info of each telemetry and science run
            print('Telemetry scan data: ')
            for il in range(len(telBegin)):
                print('Scan #' + str(il))
                print(telBegin[il], telEnd[il])

            print('Science scan data: ')
            for il in range(len(sciBegin)):
                print('Scan #' + str(il))
                print(sciBegin[il], sciEnd[il])

        #Transform the data to ndarray(np.array)
        amp = np.array(amp, dtype = object)
        tempSipm = np.array(tempSipm)
        tempAdc = np.array(tempAdc)
        vMon = np.array(vMon)
        iMon = np.array(iMon)
        #iMon = iMon / 2.0 #Adjustment for new analog circuit, from 2k resistor to 1k resistor
        bias = np.array(bias)
        uscount = np.array(uscount)
        utc = np.array(utc)
        uscountEvt = np.array(uscountEvt, dtype = object)
        if not rateStyle == '':
            timeCorrect = np.array(timeCorrect)
            timeCorrect = timeCorrect[(timeCorrect > 0) * (timeCorrect < 20 * np.std(timeCorrect))]
        if newProgramme:
            effectiveCount = np.array(effectiveCount)
            missingCount = np.array(missingCount)
        if isBin:
            sciNum = np.array(sciNum, dtype = object)
            telNum = np.array(telNum)

        #Time cut
        if isBin:
            if extractFilename(filename) in parameters.cutFileRef:
                timeCut = parameters.cutFileRef[extractFilename(filename)]
        if timeCut is not None:
            if isCi == 2:
                for isc in range(len(uscount)):
                    q1 = None
                    if isinstance(timeCut, list):
                        q1 = (np.array(uscount[isc]) >= timeCut[0]) * (np.array(uscount[isc]) <= timeCut[1])
                    else:
                        q1 = np.array(uscount[isc]) >= timeCut
                    uscount[isc] = np.array(uscount[isc])[q1]
                    utc[isc] = np.array(utc[isc])[q1]
                    if isBin:
                        telNum[isc] = np.array(telNum[isc])[q1]
                    for ich in range(4):
                        tempSipm[ich, isc] = np.array(tempSipm[ich, isc])[q1]
                        tempAdc[ich, isc] = np.array(tempAdc[ich, isc])[q1]
                        vMon[ich, isc] = np.array(vMon[ich, isc])[q1]
                        iMon[ich, isc] = np.array(iMon[ich, isc])[q1]
                        bias[ich, isc] = np.array(bias[ich, isc])[q1]
                        q2 = None
                        if isinstance(timeCut, list):
                            q2 = (np.array(uscountEvt[ich, isc]) >= timeCut[0]) * (np.array(uscountEvt[ich, isc]) <= timeCut[1])
                        else:
                            q2 = np.array(uscountEvt[ich, isc]) >= timeCut
                        uscountEvt[ich, isc] = np.array(uscountEvt[ich, isc])[q2]
                        amp[ich, isc] = np.array(amp[ich, isc])[q2]
            else:
                q1 = None
                if isinstance(timeCut, list):
                    q1 = (uscount >= timeCut[0]) * (uscount <= timeCut[1])
                else:
                    q1 = (uscount >= timeCut)
                uscount = uscount[q1]
                utc = utc[q1]
                tempSipm = tempSipm[:, q1]
                tempAdc = tempAdc[:, q1]
                vMon = vMon[:, q1]
                iMon = iMon[:, q1]
                bias = bias[:, q1]
                if isBin:
                    telNum = telNum[q1]
                for ich in range(4):
                    q2 = None
                    if isinstance(timeCut, list):
                        q2 = (np.array(uscountEvt[ich]) >= timeCut[0]) * (np.array(uscountEvt[ich]) <= timeCut[1])
                    else:
                        q2 = np.array(uscountEvt[ich]) >= timeCut
                    uscountEvt[ich] = np.array(uscountEvt[ich])[q2]
                    amp[ich] = np.array(amp[ich])[q2]
                    if isBin:
                        sciNum[ich] = np.array(sciNum[ich])[q2]
        
        #Energy cut
        if energyCut is not None:
            if isCi == 2:
                for isc in range(len(uscount)):
                    curcorr = tempBiasCorrection(tempSipm[:, isc], bias[:, isc])[0]
                    curenergy = None
                    for ich in range(4):
                        curenergy = ecCorrection(np.array(amp[ich][isc]) * curcorr[ich], True, ich)
                        if isinstance(energyCut, list):
                            q1 = (np.array(curenergy) >= energyCut[0]) * (np.array(curenergy) <= energyCut[1])
                        else:
                            q1 = np.array(curenergy) >= energyCut
                        uscountEvt[ich][isc] = np.array(uscountEvt[ich][isc])[q1]
                        amp[ich][isc] = np.array(amp[ich][isc])[q1]
                        if isBin:
                            sciNum[ich][isc] = np.array(sciNum[ich][isc])[q1]
            else:
                curcorr = tempBiasCorrection(tempSipm, bias)[0]
                curenergy = None
                for ich in range(4):
                    curenergy = ecCorrection(np.array(amp[ich]) * curcorr[ich], True, ich)
                    if isinstance(energyCut, list):
                        q1 = (np.array(curenergy) >= energyCut[0]) * (np.array(curenergy) <= energyCut[1])
                    else:
                        q1 = np.array(curenergy) >= energyCut
                    uscountEvt[ich] = np.array(uscountEvt[ich])[q1]
                    amp[ich] = np.array(amp[ich])[q1]
                    if isBin:
                        sciNum[ich] = np.array(sciNum[ich])[q1]

        #Check for error temperature data
        if parameters.tempCheck:
            if isCi == 2:
                for isc in range(len(uscount)):
                    qTemp = np.array(tempSipm[0][isc]) > -np.inf
                    for ich in range(4):
                        curtempSiPM = list(tempSipm[ich][isc])
                        curtempSiPM.append(tempSipm[ich][isc][-1])
                        curtempSiPM = np.array(curtempSiPM)
                        qTemp *= (abs(curtempSiPM[1:] - curtempSiPM[:-1]) <= np.sqrt(np.std(tempSipm[ich][isc]) ** 2 + parameters.tempErrSys ** 2))
                    uscount[isc] = np.array(uscount[isc])[qTemp]
                    utc[isc] = np.array(utc[isc])[qTemp]
                    if isBin:
                        telNum[isc] = np.array(telNum[isc])[qTemp]
                    for ich in range(4):
                        tempSipm[ich, isc] = np.array(tempSipm[ich, isc])[qTemp]
                        tempAdc[ich, isc] = np.array(tempAdc[ich, isc])[qTemp]
                        vMon[ich, isc] = np.array(vMon[ich, isc])[qTemp]
                        iMon[ich, isc] = np.array(iMon[ich, isc])[qTemp]
                        bias[ich, isc] = np.array(bias[ich, isc])[qTemp]
            else:
                qTemp = np.array(tempSipm[0]) > -np.inf
                for ich in range(4):
                    curtempSiPM = list(tempSipm[ich])
                    curtempSiPM.append(tempSipm[ich][-1])
                    curtempSiPM = np.array(curtempSiPM)
                    qTemp *= (abs(curtempSiPM[1:] - curtempSiPM[:-1]) <= np.sqrt(np.std(tempSipm[ich]) ** 2 + parameters.tempErrSys ** 2))
                uscount = uscount[qTemp]
                utc = utc[qTemp]
                tempSipm = tempSipm[:, qTemp]
                tempAdc = tempAdc[:, qTemp]
                vMon = vMon[:, qTemp]
                iMon = iMon[:, qTemp]
                bias = bias[:, qTemp]
                if isBin:
                    telNum = telNum[qTemp]

    #------------------------------------------------------------Data output------------------------------------------------------------
    print('Data readout of \"' + extractFilename(filename) + '\" complete')

    #Create output dictionaries
    sciExtracted = {
        'amp':                     amp, 
        'timestampEvt':       uscountEvt, 
        'timeCorrect':          timeCorrect, 
        'sciNum':                 sciNum, 
    }
    if newProgramme:
        sciExtracted.update({
            'effectiveCount':      effectiveCount, 
            'missingCount':        missingCount, 
        })
    if isCi > 0:
        sciExtracted.update({
            'ampCI':                       np.array(ampCI, dtype = object), 
            'timestampEvtCI':         np.array(uscountEvtCI, dtype = object), 
        })
        if newProgramme:
            sciExtracted.update({
                'effectiveCountCI':       np.array(effectiveCountCI), 
                'missingCountCI':          np.array(missingCountCI)
            })
    
    telExtracted = {
        'tempSipm':              tempSipm, 
        'tempAdc':               tempAdc, 
        'vMon':                     vMon, 
        'iMon':                      iMon, 
        'bias':                       bias, 
        'timestamp':             uscount, 
        'utc':                        utc, 
        'telNum':                  telNum, 
    }

    #Output
    if isScan:
        scanExtracted = {
            'vSet':         np.array(vSet), 
            'vScan':       np.array(vScan), 
            'iScan':        np.array(iScan), 
        }
        return sciExtracted, telExtracted, scanExtracted
    else:
        return sciExtracted, telExtracted

def HPGeDataReadout(filename):

    """
    Function for reading out single HPGe raw outout file

    Parameters
    ----------
    filename : str,
        name of the input file, currently supporting only text(.txt) files\n

    Returns
    ----------
    time : float,
        live time extracted from data file\n
    cts : array,
        spectrum extracted from data file\n
    """

    with open(filename) as f:
        #Remove \r\n using the following line
        lines = [line.rstrip() for line in f]

    cts = []
    for line in lines:
        cts.append(int(line))
    cts = np.array(cts)
    time = float(cts[0])
    cts[0], cts[1] = 0, 0

    return time, cts

def splitData(filename, sciData, telData, singlech = False, channel = -1, splitByRun = True, timeCut = None, energyCut = None):

    """
    NOTE: FOR BINARY MULTIPLE RUNS(TEMP-BIAS SCAN) DATA ONLY
    ---------

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
                    'timestampEvt':         array of arrays, extracted event uscount values, which corresponds to all events recorded in 'amp',
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event,
                    'timeCorrect':           array of arrays, correct live time, calculated with event uscount values, if rate correction is reqiured, 
                }
    telData : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telData = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'tempAdc':                array of arrays, extracted ADC temperature, 
                    'vMon':                     array of arrays, extracted monitored voltage, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values(bias = vMon - iMon * R), 
                    'timestamp':             array-like or array of arrays, extracted uscount values, 
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
                    'sciNum':                   array of arrays, split science run number of each recorded event,
                    'timeCorrect':           array of arrays, same as 'timeCorrect' in original dictioary sciData, if rate correction is reqiured, 
                }
    telDataFinal : dict
        split telemetry data, in the form of a dictionary as\n
        telDataFinal = {
                    'tempSipm':              array of arrays, split SiPM temperature, in similar form as 'amp' in sciDataFinal, 
                    'tempAdc':                array of arrays, split ADC temperature, in the same form as 'tempSipm', 
                    'vMon':                     array of arrays, split monitored voltage, in the same form as 'tempSipm', 
                    'iMon':                      array of arrays, split monitored current, in the same form as 'tempSipm', 
                    'bias':                       array of arrays, split bias, in the same form as 'tempSipm', 
                    'timestamp':             array-like or array of arrays, split uscount values, in the form of `array(array, len = NumOfRuns)`, 
                    'utc':                        array-like or array of arrays, split UTC readout value, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, split telemetry run number of each temetry data, in the same form as 'timestamp', 
                }

    Raises
    ----------
    Execption
        when the data is single-channeled and given channel number is out of bound [0-3]
    """

    if singlech:
        if not isChannel(channel):
            raise Exception('splitScans: channel number out of bound[0-3]')
    
    #Read each data from extracted data list
    amp = sciData['amp']
    uscountEvt = sciData['timestampEvt']
    sciNum = sciData['sciNum']
    tempSipm = telData['tempSipm']
    tempAdc = telData['tempAdc']
    vMon = telData['vMon']
    iMon = telData['iMon']
    bias = telData['bias']
    uscount = telData['timestamp']
    utc = telData['utc']
    telNum = telData['telNum']

    #Construct new data lists
    ampFinal = []
    uscountEvtFinal = []
    sciNumFinal = []
    tempSipmFinal = []
    tempAdcFinal = []
    vMonFinal = []
    iMonFinal = []
    biasFinal = []
    uscountFinal = []
    utcFinal = []
    telNumFinal = []
    for ich in range(4):
        ampFinal.append([])
        uscountEvtFinal.append([])
        sciNumFinal.append([])
        tempSipmFinal.append([])
        tempAdcFinal.append([])
        biasFinal.append([])
        vMonFinal.append([])
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
                if filename in parameters.unusedTelRef and curTelNum in parameters.unusedTelRef[filename]:
                    curTelNum += 1
                    continue
                qTel = np.array(telNum) == curTelNum
                uscountFinal.append(list(np.array(uscount)[qTel]))
                utcFinal.append(list(np.array(utc)[qTel]))
                telNumFinal.append(list(np.array(telNum)[qTel]))
                for ich in range(4):
                    tempSipmFinal[ich].append(list(np.array(tempSipm[ich])[qTel]))
                    tempAdcFinal[ich].append(list(np.array(tempAdc[ich])[qTel]))
                    vMonFinal[ich].append(list(np.array(vMon[ich])[qTel]))
                    iMonFinal[ich].append(list(np.array(iMon[ich])[qTel]))
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
                if filename in parameters.unusedSciRef and curSciNum in parameters.unusedSciRef[filename]:
                    curSciNum += 1
                    continue
                qSci = np.array(sciNum[ich]) == curSciNum
                ampFinal[ich].append(list(np.array(amp[ich])[qSci]))
                uscountEvtFinal[ich].append(list(np.array(uscountEvt[ich])[qSci]))
                sciNumFinal[ich].append(list(np.array(sciNum[ich])[qSci]))
                curSciNum += 1
        
        #Time cut
        if timeCut is not None:
            #Telemetry data cut
            for itel in range(len(uscountFinal)):
                q1 = None
                if isinstance(timeCut, list):
                    q1 = (np.array(uscountFinal[itel]) >= timeCut[0]) * (np.array(uscountFinal[itel]) <= timeCut[1])
                else:
                    q1 = np.array(uscountFinal[itel]) >= timeCut
                    uscountFinal[itel] = list(np.array(uscountFinal[itel])[q1])
                    utcFinal[itel] = list(np.array(utcFinal[itel])[q1])
                    telNumFinal[itel] = list(np.array(telNumFinal[itel])[q1])
                    for ich in range(4):
                        tempSipmFinal[ich][itel] = list(np.array(tempSipmFinal[ich][itel])[q1])
                        tempAdcFinal[ich][itel] = list(np.array(tempAdcFinal[ich][itel])[q1])
                        vMonFinal[ich][itel] = list(np.array(vMonFinal[ich][itel])[q1])
                        iMonFinal[ich][itel] = list(np.array(iMonFinal[ich][itel])[q1])
                        biasFinal[ich][itel] = list(np.array(biasFinal[ich][itel])[q1])
            
            #Science data cut
            for ich in range(4):
                for isci in range(len(uscountEvtFinal[ich])):
                    q2 = None
                    if isinstance(timeCut, list):
                        q2 = (np.array(uscountEvtFinal[ich][isci]) >= timeCut[0]) * (np.array(uscountEvtFinal[ich][isci]) <= timeCut[1])
                    else:
                        q2 = np.array(uscountEvtFinal[ich][isci]) >= timeCut
                    ampFinal[ich][isci] = list(np.array(ampFinal[ich][isci])[q2])
                    uscountEvtFinal[ich][isci] = list(np.array(uscountEvtFinal[ich][isci])[q2])
                    sciNumFinal[ich][isci] = list(np.array(sciNumFinal[ich][isci])[q2])
        
        #Energy cut
        if energyCut is not None:
            tempCorr = np.array(tempSipmFinal, dtype = object)
            biasCorr = np.array(biasFinal, dtype = object)
            for ich in range(4):
                if len(uscountFinal) != len(uscountEvtFinal[ich]):
                    print('splitData: Unable to conduct energy cut due to non-corresponding number of science and telemetry runs. Cut session aborted')
                    break
                for isci in range(len(uscountEvtFinal[ich])):
                    curcorr = tempBiasCorrection(tempCorr[:, isci], biasCorr[:, isci])[0]
                    curenergy = None
                    for ich in range(4):
                        curenergy = ecCorrection(np.array(ampFinal[ich][isci]) * curcorr[ich], True, ich)
                        if isinstance(energyCut, list):
                            q1 = (np.array(curenergy) >= energyCut[0]) * (np.array(curenergy) <= energyCut[1])
                        else:
                            q1 = np.array(curenergy) >= energyCut
                        uscountEvtFinal[ich][isci] = list(np.array(uscountEvtFinal[ich][isci])[q1])
                        ampFinal[ich][isci] = list(np.array(ampFinal[ich][isci])[q1])
                        sciNumFinal[ich][isci] = list(np.array(sciNumFinal[ich][isci])[q1])
    
    #--------------------------------------------------Split data according to specified bias values--------------------------------------------------
    else:
        #Get bias reference values and cut boundary for the data cut
        refVaules = parameters.biasRef['bias_set']
        bound = parameters.biasRef['bound']
        if filename in parameters.unusedRef:
            #Skip unused values in some files
            refVaules = parameters.unusedRef[filename]
        print('splitScans: Splitting data according to set bias ' + str(refVaules) + ' with bounds ' + str(bound))
        for biasSet in refVaules:
            timeBegin = -1
            timeEnd = np.max(uscount) + 1
            #Record beginning and ending time of each split runs
            if singlech:
                qBias = (np.array(bias[channel]) >= (biasSet - bound)) * (np.array(bias[channel]) <= (biasSet + bound))
                if timeBegin < np.min(np.array(uscount)[qBias]):
                    timeBegin = np.min(np.array(uscount)[qBias])
                if timeEnd > np.max(np.array(uscount)[qBias]):
                    timeEnd = np.max(np.array(uscount)[qBias])
            else:
                for ich in range(4):
                    qBias = (np.array(bias[ich]) >= (biasSet - bound)) * (np.array(bias[ich]) <= (biasSet + bound))
                    if timeBegin < np.min(np.array(uscount)[qBias]):
                        timeBegin = np.min(np.array(uscount)[qBias])
                    if timeEnd > np.max(np.array(uscount)[qBias]):
                        timeEnd = np.max(np.array(uscount)[qBias])
            #Split the data according to the beginning and ending times
            qTel = (np.array(uscount) >= timeBegin) * (np.array(uscount) <= timeEnd)
            uscountFinal.append(list(np.array(uscount)[qTel]))
            utcFinal.append(list(np.array(utc)[qTel]))
            telNumFinal.append(list(np.array(telNum)[qTel]))
            for ich in range(4):
                tempSipmFinal[ich].append(list(np.array(tempSipm[ich])[qTel]))
                tempAdcFinal[ich].append(list(np.array(tempAdc[ich])[qTel]))
                vMonFinal[ich].append(list(np.array(vMon[ich])[qTel]))
                iMonFinal[ich].append(list(np.array(iMon[ich])[qTel]))
                biasFinal[ich].append(list(np.array(bias[ich])[qTel]))
                #Split science data
                qSci = (np.array(uscountEvt[ich]) >= timeBegin) * (np.array(uscountEvt[ich]) <= timeEnd)
                ampFinal[ich].append(list(np.array(amp[ich])[qSci]))
                uscountEvtFinal[ich].append(list(np.array(uscountEvt[ich])[qSci]))
                sciNumFinal[ich].append(list(np.array(sciNum[ich])[qSci]))
    
    #Construct new data dictionaries with split data
    sciDataFinal = {
        'amp':                      np.array(ampFinal, dtype = object), 
        'timestampEvt':        np.array(uscountEvtFinal, dtype = object), 
        'sciNum':                  np.array(sciNumFinal, dtype = object), 
    }
    if 'timeCorrect' in sciData:
        sciDataFinal.update({
            'timeCorrect':           sciData['timeCorrect'], 
        })
    telDataFinal = {
        'tempSipm':             np.array(tempSipmFinal, dtype = object), 
        'tempAdc':               np.array(tempAdcFinal, dtype = object), 
        'vMon':                     np.array(vMonFinal, dtype = object), 
        'iMon':                      np.array(iMonFinal, dtype = object), 
        'bias':                       np.array(biasFinal, dtype = object), 
        'timestamp':             np.array(uscountFinal, dtype = object), 
        'utc':                        np.array(utcFinal, dtype = object), 
        'telNum':                  np.array(telNumFinal, dtype = object), 
    }
    return sciDataFinal, telDataFinal

def deleteEmptyScan(sciData, telData, scanRange = [], rateStyle = '', newProgramme = False):

    """
    Auxiliary function to delete empty sublists for ranged multiple scan data

    Parameters
    ----------
    sciData : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciData = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event uscount values, which corresponds to all events recorded in amp,
                    'timeCorrect':           array of arrays, correct live time, calculated with event uscount values, if rate correction is reqiured, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 
                    'effectiveCount':       array of arrays, extracted effective count values, for new programme files only, 
                    'missingCount':          array of arrays, extracted missing count values, for new programme files only, 

                    for data with charge injection (CI) part:
                    'ampCI':                     array of arrays, extracted amplitude data of CI part, 
                    'timestampEvtCI':       array of arrays, extracted event uscount values of CI part, which corresponds to all events recorded in ampCI, 
                    'effectiveCountCI':     array of arrays, extracted effective count values of CI part, 
                    'missingCountCI':        array of arrays, extracted missing count values of CI part, 
                }
    telData : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telData = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'tempAdc':                array of arrays, extracted ADC temperature, 
                    'vMon':                     array of arrays, extracted monitored voltage, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values(bias = vMon - iMon * R), 
                    'timestamp':             array-like or array of arrays, extracted uscount values, 
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
                    'timestampEvt':         array of arrays, extracted event uscount values, which corresponds to all events recorded in amp,
                    'timeCorrect':           array of arrays, correct live time, calculated with event uscount values, if rate correction is reqiured, 
                    'effectiveCount':       array of arrays, extracted effective count values, for new programme files only, 
                    'missingCount':          array of arrays, extracted missing count values, for new programme files only, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event,

                    for data with charge injection (CI) part:
                    'ampCI':                     array of arrays, extracted amplitude data of CI part, 
                    'timestampEvtCI':       array of arrays, extracted event uscount values of CI part, which corresponds to all events recorded in ampCI, 
                    'effectiveCountCI':     array of arrays, extracted effective count values of CI part, for new programme files only, 
                    'missingCountCI':        array of arrays, extracted missing count values of CI part, for new programme files only, 
                }
    telDataCleaned : dict
        telemetry data with empty runs removed, in the form of a dictionary as\n
        telDataCleaned = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'tempAdc':                array of arrays, extracted ADC temperature, 
                    'vMon':                     array of arrays, extracted monitored voltage, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values(bias = vMon - iMon * R), 
                    'timestamp':             array of arrays, extracted uscount values, 
                    'utc':                        array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }
    """

    #Read each data from extracted data list
    amp = sciData['amp']
    uscountEvt = sciData['timestampEvt']
    sciNum = sciData['sciNum']
    timeCorrect = []
    if not rateStyle == '':
        timeCorrect = sciData['timeCorrect']
    tempSipm = telData['tempSipm']
    tempAdc = telData['tempAdc']
    vMon = telData['vMon']
    iMon = telData['iMon']
    bias = telData['bias']
    uscount = telData['timestamp']
    utc = telData['utc']
    telNum = telData['telNum']
    effectiveCount = []
    missingCount = []
    ampCI = []
    uscountEvtCI = []
    effectiveCountCI = []
    missingCountCI = []
    if 'effectiveCount' in sciData:
        effectiveCount = sciData['effectiveCount']
    if 'missingCount' in sciData:
        missingCount = sciData['missingCount']
    if 'ampCI' in sciData:
        ampCI = sciData['ampCI']
    if 'timestampEvtCI' in sciData:
        uscountEvtCI = sciData['timestampEvtCI']
    if 'effectiveCountCI' in sciData:
        effectiveCountCI = sciData['effectiveCountCI']
    if 'missingCountCI' in sciData:
        missingCountCI = sciData['missingCountCI']

    if len(scanRange) == 0:
        lower = 0
        upper = len(uscount) - 1
    else:
        lower = scanRange[0] - 1
        upper = scanRange[1] - 1
    nScan = len(uscount)

    #Convert all data to list of lists in order to delete empty scans
    amp = list(amp)
    uscountEvt = list(uscountEvt)
    sciNum = list(sciNum)
    tempSipm = list(tempSipm)
    tempAdc = list(tempAdc)
    vMon = list(vMon)
    iMon = list(iMon)
    bias = list(bias)
    uscount = list(uscount)
    utc = list(utc)
    telNum = list(telNum)
    if not rateStyle == '':
        timeCorrect = list(timeCorrect)
    ampCI = list(ampCI)
    uscountEvtCI = list(uscountEvtCI)
    if newProgramme:
        effectiveCount = list(effectiveCount)
        effectiveCountCI = list(effectiveCountCI)
        missingCount = list(missingCount)
        missingCountCI = list(missingCountCI)
    for isc in range(nScan):
        uscount[isc] = list(uscount[isc])
        utc[isc] = list(utc[isc])
        if not rateStyle == '':
            timeCorrect[isc] = list(timeCorrect[isc])
        if len(telNum) > 0:
            telNum[isc] = list(telNum[isc])
        if newProgramme:
            if len(effectiveCount) > 0:
                effectiveCount[isc] = list(effectiveCount[isc])
                missingCount[isc] = list(missingCount[isc])
            if len(effectiveCountCI) > 0:
                effectiveCountCI[isc] = list(effectiveCountCI[isc])
                missingCountCI[isc] = list(missingCountCI[isc])
    for ich in range(4):
        amp[ich] = list(amp[ich])
        tempSipm[ich] = list(tempSipm[ich])
        tempAdc[ich] = list(tempAdc[ich])
        vMon[ich] = list(vMon[ich])
        iMon[ich] = list(iMon[ich])
        bias[ich] = list(bias[ich])
        uscountEvt[ich] = list(uscountEvt[ich])
        if len(ampCI) > 0:
            ampCI[ich] = list(ampCI[ich])
            uscountEvtCI[ich] = list(uscountEvtCI[ich])
        if len(sciNum) > 0:
            sciNum[ich] = list(sciNum[ich])
        for isc in range(nScan):
            amp[ich][isc] = list(amp[ich][isc])
            tempSipm[ich][isc] = list(tempSipm[ich][isc])
            tempAdc[ich][isc] = list(tempAdc[ich][isc])
            vMon[ich][isc] = list(vMon[ich][isc])
            iMon[ich][isc] = list(iMon[ich][isc])
            bias[ich][isc] = list(bias[ich][isc])
            uscountEvt[ich][isc] = list(uscountEvt[ich][isc])
            if len(ampCI) > 0 and len(ampCI[ich]) > 0:
                ampCI[ich][isc] = list(ampCI[ich][isc])
                uscountEvtCI[ich][isc] = list(uscountEvtCI[ich][isc])
            if len(sciNum[ich]) > 0:
                sciNum[ich][isc] = list(sciNum[ich][isc])

    #Delete scans before first required scan
    for isc in range(lower):
        del uscount[0]
        del utc[0]
        if not rateStyle == '':
            del timeCorrect[0]
        if len(telNum) > 0:
            del telNum[0]
        if newProgramme:
            if len(effectiveCount) > 0:
                del effectiveCount[0], missingCount[0]
            if len(effectiveCountCI) > 0:
                del effectiveCountCI[0], missingCountCI[0]
        for ich in range(4):
            del amp[ich][0], tempSipm[ich][0], tempAdc[ich][0], vMon[ich][0], iMon[ich][0], bias[ich][0], uscountEvt[ich][0]
            if len(ampCI) > 0 and len(ampCI[ich]) > 0:
                del ampCI[ich][0], uscountEvtCI[ich][0]
            if len(sciNum) > 0 and len(sciNum[ich]) > 0:
                del sciNum[ich][0]

    #Delete scans after last required scan
    for isc in range(nScan - upper - 1):
        del uscount[-1]
        del utc[-1]
        if not rateStyle == '':
            del timeCorrect[-1]
        if len(telNum) > 0:
            del telNum[-1]
        if newProgramme:
            if len(effectiveCount) > 0:
                del effectiveCount[-1], missingCount[-1]
            if len(effectiveCountCI) > 0:
                del effectiveCountCI[-1], missingCountCI[-1]
        for ich in range(4):
            del amp[ich][-1], tempSipm[ich][-1], tempAdc[ich][-1], vMon[ich][-1], iMon[ich][-1], bias[ich][-1], uscountEvt[ich][-1]
            if len(ampCI) > 0 and len(ampCI[ich]) > 0:
                del ampCI[ich][-1], uscountEvtCI[ich][-1]
            if len(sciNum) > 0 and len(sciNum[ich]) > 0:
                del sciNum[ich][-1]
    
    #Delete other empty scans/runs
    currPos = 0
    while currPos < len(uscount):
        if len(uscount[currPos]) == 0:
            del uscount[currPos], utc[currPos]
            if not rateStyle == '':
                del timeCorrect[currPos]
            if len(telNum) > 0:
                del telNum[currPos]
            if newProgramme:
                if len(effectiveCount) > 0:
                    del effectiveCount[currPos], missingCount[currPos]
                if len(effectiveCountCI) > 0:
                    del effectiveCountCI[currPos], missingCountCI[currPos]
            for ich in range(4):
                del amp[ich][currPos], uscountEvt[ich][currPos], tempSipm[ich][currPos], tempAdc[ich][currPos], iMon[ich][currPos], vMon[ich][currPos], \
                    bias[ich][currPos]
                if len(ampCI) > 0 and len(ampCI[ich]) > 0:
                    del ampCI[ich][currPos], uscountEvtCI[ich][currPos]
                if len(sciNum) > 0 and len(sciNum[ich]) > 0:
                    del sciNum[ich][currPos]
        else:
            currPos += 1
    
    #Convert data back to ndarray form and create output dictionaries
    sciDataCleaned = {
        'amp':                      np.array(amp, dtype = object), 
        'timestampEvt':        np.array(uscountEvt, dtype = object), 
        'sciNum':                  np.array(sciNum, dtype = object), 
    }
    if not rateStyle == '':
        sciDataCleaned.update({
            'timeCorrect':          np.array(timeCorrect, dtype = object), 
        })
    if len(effectiveCount) > 0:
        sciDataCleaned.update({
            'effectiveCount':           np.array(effectiveCount), 
            'missingCount':              np.array(missingCount), 
        })
    if len(ampCI) > 0:
        sciDataCleaned.update({
            'ampCI':                        np.array(ampCI, dtype = object), 
            'timestampEvtCI':          np.array(uscountEvtCI, dtype = object), 
        })
    if len(effectiveCountCI) > 0:
        sciDataCleaned.update({
            'effectiveCountCI':           np.array(effectiveCountCI, dtype = object), 
            'missingCountCI':              np.array(missingCountCI, dtype = object), 
        })
    telDataCleaned = {
        'tempSipm':              np.array(tempSipm, dtype = object), 
        'tempAdc':                np.array(tempAdc, dtype = object), 
        'vMon':                     np.array(vMon, dtype = object), 
        'iMon':                      np.array(iMon, dtype = object), 
        'bias':                       np.array(bias, dtype = object), 
        'timestamp':             np.array(uscount, dtype = object), 
        'utc':                        np.array(utc, dtype = object), 
        'telNum':                  np.array(telNum, dtype = object), 
    }

    #Thanks to ndarray, ALL THESE matter is needed to delete a SINGLE sublist!
    return sciDataCleaned, telDataCleaned

def fileOutput(file, sciData, telData, scanData = {}, isCi = 0, isScan = False, rateStyle = '', newProgramme = False):

    """
    NOTE: CURRENTLY USED FOR NON-MUMTIPLE-SCANS DATA ONLY
    ----------

    Function for writing designated readout data to files(numpy form binary file .npy)

    Parameters
    ----------
    file : str,
        name of the raw data file\n
    sciData : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciData = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event uscount values, which corresponds to all events recorded in amp,
                    'timeCorrect':           array of arrays, correct live time, calculated with event uscount values, if rate correction is reqiured, 
                    'effectiveCount':       array of arrays, extracted effective count values, for new programme files only, 
                    'missingCount':          array of arrays, extracted missing count values, for new programme files only, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 

                    for data with charge injection (CI) part:
                    'ampCI':                     array of arrays, extracted amplitude data of CI part, 
                    'timestampEvtCI':       array of arrays, extracted event uscount values of CI part, which corresponds to all events recorded in ampCI, 
                    'effectiveCountCI':     array of arrays, extracted effective count values of CI part, for new programme files only, 
                    'missingCountCI':        array of arrays, extracted missing count values of CI part, for new programme files only, 
                }
    telData : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telData = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'tempAdc':                array of arrays, extracted ADC temperature, 
                    'vMon':                     array of arrays, extracted monitored voltage, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values(bias = vMon - iMon * R), 
                    'timestamp':             array-like or array of arrays, extracted uscount values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', 
                }
    scanData : dict, optional
        I-V scan data extracted from data packs, in the form of a dictionary as\n
        scanData = {
                    'vSet':                      array-like, set voltage, 
                    'vScan':                    array of arrays, monitored I-V scan voltages, in the form of array(array, array, array, array), 
                    'iScan':                     array of arrays, monitored I-V scan currents, in the form of array(array, array, array, array), 
                }
    isCi : int, optional
        indicates whether the input file has CI part, with 0 for no CI, 1 for with CI\n
    isScan : boolean, optional
        indicates whether the input file has I-V scan part\n
    rateStyle : str, optional
        the style of count rate correction, '' for no corraction, 's' for correction with small data packs(512byte), 'l' for calculation with large data \
            packs(4096byte)\n
    newProgramme : boolean, optional
        indicates whether the data is produced with new hardware programme(6th ver.)\n

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
    Science run number: `sciNum_filename.npy`\n
    SiPM temperature: `tempSiPM_filename.npy`\n
    ADC temperature: `tempADC_filename.npy`\n
    Monitored voltage: `vMon_filename.npy`\n
    Monitored current:  `iMon_filename.npy`\n
    Bias: `bias_filename.npy`\n
    Timestamp: `timestamp_filename.npy`\n
    UTC: `utc_filename.npy`\n
    Telemetry run number: `telNum_filename.npy`\n
    Correct live time: `livetime_filename.npy`\n
    Effective count: `eff_filename.npy`\n
    Missing count: `miss_filename.npy`\n
    ----
    CI data(output when `isCi != 0`): \n
    CI ADC amplitude: `ci_amp_filename.npy`\n
    CI event uscount: `ci_uscevt_filename.npy`\n
    CI effective count: `ci_eff_filename.npy`\n
    CI missing count: `ci_miss_filename.npy`\n
    ----
    Scan data: (output when `isScan == True`)\n
    Set voltage for voltage scan: `vSet_filename.npy`\n
    Scan voltage for voltage scan: `vScan_filename.npy`\n
    Scan current for voltage scan: `iScan_filename.npy`\n
    """

    amp = sciData['amp']
    sciNum = sciData['sciNum']
    timestampEvt = sciData['timestampEvt']
    timeCorrect = []
    if not rateStyle == '':
        timeCorrect = sciData['timeCorrect']
    if newProgramme:
        effectiveCount = sciData['effectiveCount']
        missingCount = sciData['missingCount']
    tempSipm = telData['tempSipm']
    tempAdc = telData['tempAdc']
    vMon = telData['vMon']
    iMon = telData['iMon']
    bias = telData['bias']
    timestamp = telData['timestamp']
    utc = telData['utc']
    telNum = telData['telNum']
    ampCI = []
    timestampEvtCI = []
    effectiveCountCI = []
    missingCountCI = []
    if isCi > 0:
        ampCI = sciData['ampCI']
        timestampEvtCI = sciData['timestampEvtCI']
        if newProgramme:
            effectiveCountCI = sciData['effectiveCountCI']
            missingCountCI = sciData['missingCountCI']
    if isScan:
        vSet = scanData['vSet']
        vScan = scanData['vScan']
        iScan = scanData['iScan']
    
    print('fileOutput: writing output files of ' + file)
    outputSuccess = True
    errCount = 0

    #Read config file
    config = {}
    try:
        with open(parameters.outputConfig, 'r') as fin:
            config = json.load(fin)
    except:
        print('fileOutput: unable to read file output configuration from file \"' + parameters.outputConfig + '\"')
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
        if config['sciNum'] and filename.endswith('.dat'):
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

        #ADC temperature
        if config['tempAdc']:
            try:
                with open(savePath + 'tempADC_' + filename + '.npy', 'wb') as fouttempadc:
                    np.save(fouttempadc, tempAdc)
            except:
                print('fileOutput: Error writing ADC temperature file \"' + 'tempADC_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Monitored voltage
        if config['vMon']:
            try:
                with open(savePath + 'vMon_' + filename + '.npy', 'wb') as foutvMon:
                    np.save(foutvMon, vMon)
            except:
                print('fileOutput: Error writing monitored voltage file \"' + 'vMon_'+ filename + '.npy' + '\"')
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

        #Uscount
        if config['timestamp']:
            try:
                with open(savePath + 'timestamp_' + filename + '.npy', 'wb') as fouttime:
                    np.save(fouttime, timestamp)
            except:
                print('fileOutput: Error writing timestamp file \"' + 'usc_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #UTC
        if config['utc']:
            try:
                with open(savePath + 'utc_' + filename + '.npy', 'wb') as fouttime:
                    np.save(fouttime, utc)
            except:
                print('fileOutput: Error writing utc file \"' + 'utc_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Telemetry run number
        if config['telNum'] and filename.endswith('.dat'):
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
                print('fileOutput: Error writing event timestamp file \"' + 'uscevt_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Correct live time
        if config['livetime'] and rateStyle != '':
            try:
                with open(savePath + 'livetime_' + filename + '.npy', 'wb') as fouttimecorr:
                    np.save(fouttimecorr, timeCorrect)
            except:
                print('fileOutput: Error writing correct live time file \"' + 'livetime_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Effective count
        if config['effCount'] and newProgramme:
            try:
                with open(savePath + 'eff_' + filename + '.npy', 'wb') as fouteff:
                    np.save(fouteff, effectiveCount)
            except:
                print('fileOutput: Error writing effective count file \"' + 'eff_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        #Missing count
        if config['missCount'] and newProgramme:
            try:
                with open(savePath + 'miss_' + filename + '.npy', 'wb') as foutmiss:
                    np.save(foutmiss, missingCount)
            except:
                print('fileOutput: Error writing missing count file \"' + 'miss_'+ filename + '.npy' + '\"')
                outputSuccess = False
                errCount += 1

        ind = 12

        #CI data
        if isCi > 0:
            #CI ADC amplitude
            if config['ampCI']:
                try:
                    with open('ci_amp_' + filename + '.npy', 'wb') as foutampCI:
                        np.save(foutampCI, ampCI)
                except:
                    print('fileOutput: Error writing CI ADC amplitude file file \"' + 'ci_amp_'+ filename + '.npy' + '\"')
                    outputSuccess = False
                    errCount += 1
            ind += 1

            #CI event time data
            if config['timestampEvtCI']:
                try:
                    with open('ci_timeevt_' + filename + '.npy', 'wb') as fouttimeevtCI:
                        np.save(fouttimeevtCI, timestampEvtCI)
                except:
                    print('fileOutput: Error writing CI event timestamp file \"' + 'ci_timeevt_'+ filename + '.npy' + '\"')
                    outputSuccess = False
                    errCount += 1
            ind += 1

            #CI effective count
            if config['effCountCI'] and newProgramme:
                try:
                    with open('ci_eff_' + filename + '.npy', 'wb') as fouteffCI:
                        np.save(fouteffCI, effectiveCountCI)
                except:
                    print('fileOutput: Error writing CI effective count file \"' + 'ci_eff_'+ filename + '.npy' + '\"')
                    outputSuccess = False
                    errCount += 1
            ind += 1

            #CI missing count
            if config['missCountCI'] and newProgramme:
                try:
                    with open('ci_miss_' + filename + '.npy', 'wb') as foutmissCI:
                        np.save(foutmissCI, missingCountCI)
                except:
                    print('fileOutput: Error writing CI missing count file \"' + 'ci_miss_'+ filename + '.npy' + '\"')
                    outputSuccess = False
                    errCount += 1
            ind += 1

        #I-V scan
        if isScan:
            #Set voltage
            if config['vSet']:
                try:
                    with open('vSet_' + filename + '.npy', 'wb') as foutvSet:
                        np.save(foutvSet, vSet)
                except:
                    print('fileOutput: Error writing set voltage file \"' + 'vSet_'+ filename + '.npy' + '\"')
                    outputSuccess = False
                    errCount += 1
            ind += 1

            #Scan voltage
            if config['vScan']:
                try:
                    with open('vScan_' + filename + '.npy', 'wb') as foutvScan:
                        np.save(foutvScan, vScan)
                except:
                    print('fileOutput: Error writing scan voltage file \"' + 'vScan_'+ filename + '.npy' + '\"')
                    outputSuccess = False
                    errCount += 1
            ind += 1

            #Scan current
            if config['iScan']:
                try:
                    with open('iScan_' + filename + '.npy', 'wb') as foutiScan:
                        np.save(foutiScan, iScan)
                except:
                    print('fileOutput: Error writing scan current file \"' + 'iScan_'+ filename + '.npy' + '\"')
                    outputSuccess = False
                    errCount += 1

        print('File output complete with ' + str(errCount) + ' error(s)')

    except:
        print('fileOutput: error occurred while reading config file, abandoning file output process')
        outputSuccess = False

    return outputSuccess

def importData(filename, importPath, isCi = 0, isScan = False):

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
    isScan : boolean, optional
        indicates whether the input file has I-V scan part\n

    Returns
    ----------
    sciExtracted : dict
        science data extracted from science data packs, in the form of a dictionary as\n
        sciExtracted = {
                    'amp':                       array of arrays, extracted amplitude data,
                    'timestampEvt':         array of arrays, extracted event uscount values, which corresponds to all events recorded in amp,
                    'timeCorrect':           array of arrays, correct live time, calculated with event uscount values, if rate correction is reqiured, 
                    'effectiveCount':       array of arrays, extracted effective count values, for new programme files only, 
                    'missingCount':          array of arrays, extracted missing count values, for new programme files only, 
                    'sciNum':                   array of arrays, corresponding science run number of each recorded event, for binary files only, 

                    for data with charge injection (CI) part:
                    'ampCI':                     array of arrays, extracted amplitude data of CI part, 
                    'timestampEvtCI':       array of arrays, extracted event uscount values of CI part, which corresponds to all events recorded in ampCI, 
                    'effectiveCountCI':     array of arrays, extracted effective count values of CI part, for new programme files only, 
                    'missingCountCI':        array of arrays, extracted missing count values of CI part, for new programme files only, 
                }
        All data are in the form of `array(array, array, array, array)` for single scan(including multiple runs), or `array(array(array, len = NumOfScans), \
            array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans))` for multiple scans
    telExtracted : dict
        telemetry data extracted from telemetry data packs, in the form of a dictionary as\n
        telExtracted = {
                    'tempSipm':              array of arrays, extracted SiPM temperature, 
                    'tempAdc':                array of arrays, extracted ADC temperature, 
                    'vMon':                     array of arrays, extracted monitored voltage, 
                    'iMon':                      array of arrays, extracted monitored current, 
                    'bias':                       array of arrays, extracted bias values(bias = vMon - iMon * R), 
                    'timestamp':             array-like or array of arrays, extracted uscount values, 
                    'utc':                        array-like or array of arrays, extracted UTC values, in the same form as 'timestamp', 
                    'telNum':                  array-like or array of arrays, corresponding telemetry run number of each temetry data, in the same form as 'timestamp', \
                        for binary files only, 
                }
        Channel-seperated data are in the form of `array(array, array, array, array)` for single scan(including multiple runs), or \
            `array(array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans), array(array, len = NumOfScans))` \
            for multiple scans\n
        Shared data are in the form of `array` for single scan(including multiple runs), or \
            `array(array, len = NumOfScans)` for multiple scans\n
    scanExtracted : dict, optional
        I-V scan data extracted from data packs, in the form of a dictionary as\n
        telExtracted = {
                    'vSet':                      array-like, set voltage, 
                    'vScan':                    array of arrays, monitored I-V scan voltages, in the form of array(array, array, array, array), 
                    'iScan':                     array of arrays, monitored I-V scan currents, in the form of array(array, array, array, array), 
                }
    """

    #Read config file
    config = {}
    readConfig = True
    try:
        with open(parameters.importConfig, 'r') as fin:
            config = json.load(fin)
    except:
        print('importData: unable to read file output configuration from file \"' + parameters.importConfig + '\"')
        print('WARNING: empty data sets will be returned from \'importData\'')
        readConfig = False

    importNames = []
    #Preprocessing of file names and paths
    if len(importPath) == 0:
        importPath.append(config['import_dir'])
    if (len(filename) == 0 or len(filename) == 1 and filename[0] == '') and config['use_own_dir']:
        print('importData: unable to import from unspecified directory with \'use_own_dir\' == True in config file \"' + parameters.importConfig + '\" and '
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
    sciNum = []
    ampCI = []
    tempSipm = []
    tempAdc = []
    vMon = []
    iMon = []
    bias = []
    timestamp = [] #uscount in telemetry data
    utc = []
    telNum = []
    timestampEvt = []
    timestampEvtCI = []
    vSet = [] #i-v scan data
    vScan = [] #i-v scan data
    iScan = [] #i-v scan data
    timeCorrect = []
    effectiveCount = []
    effectiveCountCI = []
    missingCount = []
    missingCountCI = []

    #Create output dictionaries
    sciExtracted = {}
    telExtracted = {}
    if isScan:
        scanExtracted = {}

    for ich in range(4):
        amp.append([])
        sciNum.append([])
        ampCI.append([])
        tempSipm.append([])
        tempAdc.append([])
        vMon.append([])
        iMon.append([])
        bias.append([])
        timestampEvt.append([])
        timestampEvtCI.append([])
        vScan.append([])
        iScan.append([])

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

                    #Event timestamp files
                    elif file.startswith('timestampEvt') and config['timestampEvt']:
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

                    #ADC temperature files
                    elif file.startswith('tempADC') and config['tempAdc']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curtempAdc = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curtempAdc[ich])):
                                            tempAdc[ich].append(curtempAdc[ich][isc])
                                    else:
                                        tempAdc[ich].append(curtempAdc[ich])
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Monitored voltage files
                    elif file.startswith('vMon') and config['vMon']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curvMon = np.load(fin, allow_pickle = True)
                                for ich in range(4):
                                    if isCi == 2:
                                        for isc in range(len(curvMon[ich])):
                                            tempSipm[ich].append(curvMon[ich][isc])
                                    else:
                                        vMon[ich].append(curvMon[ich])
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

                    #Correct live time files
                    elif file.startswith('livetime') and config['livetime']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curtimeCorrect = np.load(fin, allow_pickle = True)
                                if isCi == 2:
                                    for isc in range(len(curtimeCorrect)):
                                        timeCorrect.append(curtimeCorrect[isc])
                                else:
                                    timeCorrect.append(curtimeCorrect)
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Effective count files
                    elif file.startswith('eff') and config['effCount']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                cureffectiveCount = np.load(fin, allow_pickle = True)
                                if isCi == 2:
                                    for isc in range(len(cureffectiveCount)):
                                        effectiveCount.append(cureffectiveCount[isc])
                                else:
                                    effectiveCount.append(cureffectiveCount)
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #Missing count files
                    elif file.startswith('miss') and config['missCount']:
                        try:
                            with open(path + '/' + file, 'rb') as fin:
                                curmissingCount = np.load(fin, allow_pickle = True)
                                if isCi == 2:
                                    for isc in range(len(curmissingCount)):
                                        missingCount.append(curmissingCount[isc])
                                else:
                                    missingCount.append(curmissingCount)
                        except:
                            print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                            errCount += 1

                    #CI data
                    elif isCi >= 1:
                        #CI effective count files
                        if file.startswith('ci_eff') and config['effCountCI']:
                            try:
                                with open(path + '/' + file, 'rb') as fin:
                                    cureffectiveCountCI = np.load(fin, allow_pickle = True)
                                    if isCi == 2:
                                        for isc in range(len(cureffectiveCountCI)):
                                            effectiveCountCI.append(cureffectiveCountCI[isc])
                                    else:
                                        effectiveCountCI.append(cureffectiveCountCI)
                            except:
                                print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                                errCount += 1

                        #CI missing count files
                        elif file.startswith('ci_miss') and config['missCountCI']:
                            try:
                                with open(path + '/' + file, 'rb') as fin:
                                    curmissingCountCI = np.load(fin, allow_pickle = True)
                                    if isCi == 2:
                                        for isc in range(len(curmissingCountCI)):
                                            missingCountCI.append(curmissingCountCI[isc])
                                    else:
                                        missingCountCI.append(curmissingCountCI)
                            except:
                                print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                                errCount += 1

                        #CI ADC amplitude files
                        elif file.startswith('ci_amp') and config['ampCI']:
                            try:
                                with open(path + '/' + file, 'rb') as fin:
                                    curampCI = np.load(fin, allow_pickle = True)
                                    for ich in range(4):
                                        if isCi == 2:
                                            for isc in range(len(curampCI[ich])):
                                                ampCI[ich].append(curampCI[ich][isc])
                                        else:
                                            ampCI[ich].append(curampCI[ich])
                            except:
                                print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                                errCount += 1

                        #CI event timestamp files
                        elif file.startswith('ci_uscevt') and config['uscountEvtCI']:
                            try:
                                with open(path + '/' + file, 'rb') as fin:
                                    curuscountEvtCI = np.load(fin, allow_pickle = True)
                                    for ich in range(4):
                                        if isCi == 2:
                                            for isc in range(len(curuscountEvtCI[ich])):
                                                timestampEvtCI[ich].append(curuscountEvtCI[ich][isc])
                                        else:
                                            timestampEvtCI[ich].append(curuscountEvtCI[ich])
                            except:
                                print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                                errCount += 1

                    #I-V scan files
                    elif isScan:
                        #Set voltage files
                        if file.startswith('vSet_') and config['vSet']:
                            try:
                                with open(path + '/' + file, 'rb') as fin:
                                    curvSet = np.load(fin, allow_pickle = True)
                                    vSet.append(curvSet)
                            except:
                                print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                                errCount += 1

                        #Scan voltage files
                        elif file.startswith('vScan_') and config['vScan']:
                            try:
                                with open(path + '/' + file, 'rb') as fin:
                                    curvScan = np.load(fin, allow_pickle = True)
                                    vScan.append(curvScan)
                            except:
                                print('importData: error reading file \"' + path + '/' + file + '\" or file does not exist')
                                errCount += 1

                        #Scan current files
                        elif file.startswith('iScan_') and config['iScan']:
                            try:
                                with open(path + '/' + file, 'rb') as fin:
                                    curiScan = np.load(fin, allow_pickle = True)
                                    iScan.append(curiScan)
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
    if config['telNum']:
        if len(utc) == 1:
            telNum = telNum[0]
        telExtracted.update({
            'telNum':              np.array(telNum), 
        })
    if config['livetime']:
        if len(timeCorrect) == 1:
            timeCorrect = timeCorrect[0]
        sciExtracted.update({
            'timeCorrect':              np.array(timeCorrect), 
        })
    if config['effCount']:
        if len(effectiveCount) == 1:
            effectiveCount = effectiveCount[0]
        sciExtracted.update({
            'effectiveCount':              np.array(effectiveCount), 
        })
    if config['missCount']:
        if len(missingCount) == 1:
            missingCount = missingCount[0]
        sciExtracted.update({
            'missingCount':              np.array(missingCount), 
        })
    if config['effCountCI']:
        if len(effectiveCountCI) == 1:
            effectiveCountCI = effectiveCountCI[0]
        sciExtracted.update({
            'effectiveCountCI':              np.array(effectiveCountCI), 
        })
    if config['missCountCI']:
        if len(missingCountCI) == 1:
            missingCountCI = missingCountCI[0]
        sciExtracted.update({
            'missingCountCI':              np.array(missingCountCI), 
        })
    if config['vSet']:
        if len(vSet) == 1:
            vSet = vSet[0]
        scanExtracted.update({
            'vSet':              np.array(vSet), 
        })
    if config['vScan']:
        if len(vScan) == 1:
            vScan = vScan[0]
        scanExtracted.update({
            'vScan':              np.array(vScan), 
        })
    if config['iScan']:
        if len(iScan) == 1:
            iScan = iScan[0]
        scanExtracted.update({
            'iScan':              np.array(iScan), 
        })
    if config['amp']:
        for ich in range(4):
            if len(amp[ich]) == 1:
                amp[ich] = amp[ich][0]
        sciExtracted.update({
            'amp':              np.array(amp, dtype = object), 
        })
    if config['sciNum']:
        for ich in range(4):
            if len(sciNum[ich]) == 1:
                sciNum[ich] = sciNum[ich][0]
        sciExtracted.update({
            'sciNum':              np.array(sciNum, dtype = object), 
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
    if config['tempAdc']:
        for ich in range(4):
            if len(tempAdc[ich]) == 1:
                tempAdc[ich] = tempAdc[ich][0]
        telExtracted.update({
            'tempAdc':              np.array(tempAdc), 
        })
    if config['vMon']:
        for ich in range(4):
            if len(vMon[ich]) == 1:
                vMon[ich] = vMon[ich][0]
        telExtracted.update({
            'vMon':              np.array(vMon), 
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
    if config['ampCI']:
        for ich in range(4):
            if len(ampCI[ich]) == 1:
                ampCI[ich] = ampCI[ich][0]
        sciExtracted.update({
            'ampCI':              np.array(ampCI), 
        })
    if config['timestampEvtCI']:
        for ich in range(4):
            if len(timestampEvtCI[ich]) == 1:
                timestampEvtCI[ich] = timestampEvtCI[ich][0]
        sciExtracted.update({
            'timestampEvtCI':              np.array(timestampEvtCI), 
        })

    #Data output
    if isScan:
        return sciExtracted, telExtracted, scanExtracted
    else:
        return sciExtracted, telExtracted

#****************************************************************************************************************************************************************
#***********************************************************************Basic plot part***********************************************************************
#****************************************************************************************************************************************************************

def plotRawData(filename, amp, nbins, corr, time, singlech = False, channel = -1, rateStyle = '', rateAll = None, doCorr = True, utc = None, ecCorr = False, \
    binWidth = None):
    #TBD: Add additional option to convert ADC channel to energy
    """
    Function for plotting the processed, unfitted data

    Parameters
    ----------
    filename : str,
        name of raw data file\n
    amp : array of array or array-like,
        ADC amplitude data\n
    nbins : int,
        number of bins to be used in the spectrum, within [1, parameters.adcMax]\n
    corr : list,
        temperature-bias correction factors for the data\n
    time : float,
        time taken to take the spectrum, usually calculated with event uscount, used when `rateStyle == ''` (no count rate correction)\n
    singlech : boolean, optional
        indicates whether the data is single-channeled\n
    channel : int, optional
        channel number for single channel fits, in range [0-3]\n
    rateStyle : str, optional
        the style of calculating real count rate, '' for none, 's' for calculation with small data pack(512byte) live time data, 'p' for calculation with large data \
            packs(4096byte)\n
    rateAll : float or list, optional
        correct count rate of all spectrum in all 4 channels calculated when reading data, or count rates for 4 seperate channels, only used when rateStyle is \
        's' or 'p'\n
    doCorr : boolean, optional
        indicates whether the temperature-bias correction will be done, to avoid warning info in the output\n
    utc : array-like, optional
        UTC readout value, to specify time period of each run for in-orbit data quick-view\n
    ecCorr : boolean, optional
        indicates whether the EC correction will be done\n
    binWidth : int or None, optional
        bin width to be used in the spectrum, within [1, parameters.adcMax], or None if not used\n

    Returns
    ----------
    Nothing

    Raises
    ----------
    Exception
        when rateStyle given is not in the available styles, or correct count rate not specified if rate correction is required, \
        also when given channel number is out of range [0, 3] when doing single channel plots
    """

    #Save plot option
    filenameNoAppend = filename[0:filename.find(filename.split('.')[-1]) - 1]
    saveFigPath = parameters.saveFigPath + '/'
    if parameters.saveFigNoPlot:
        if not os.path.isdir(saveFigPath):
            os.makedirs(saveFigPath)

    #Dummy temp-bias correction factor
    if not doCorr:
        corr = [1.0, 1.0, 1.0, 1.0]

    #For count rate correction, correct count rate calculated with fitRateCorrect should be specified
    styleAvailable = ['', 's', 'p']
    if not rateStyle in styleAvailable:
        raise Exception('plotRawData: unknown count rate correction style')
    elif rateAll is None:
        raise Exception('plotRawData: corrected count rate not specified')

    specWidth = None
    #Single channel plots
    if singlech:
        if not isChannel(channel):
            raise Exception('plotRawData: channel number out of bound[0-3]')
        dataCorr = np.array(amp[channel]) * corr[channel]
        specRange = [0., parameters.adcMax * corr[channel]]
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech, channel)
            specRange = [0., parameters.energyPlotUpper]
        specWidth = specRange[1] - specRange[0]
        spectrum, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth)
        #Divide bin count by bin width, so that the integral is uniformed to total count
        spectrum = spectrum * nbins / specWidth
    #Multiple(all-4) channel plots
    else:
        dataCorr = []
        specRange = []
        specWidth = []
        for ich in range(4):
            amp[ich] = np.array(amp[ich])
            dataCorr.append(amp[ich] * corr[ich])
            specRange.append([0., parameters.adcMax * corr[ich]])
            specWidth.append(specRange[ich][1] - specRange[ich][0])
        dataCorr = np.array(dataCorr, dtype = object)
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech)
            specRange = [0., parameters.energyPlotUpper]
            specWidth = specRange[1] - specRange[0]
        spectrum, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth)
        for ich in range(4):
            #Divide bin count by bin width, so that the integral is uniformed to total count
            if isinstance(specWidth, list):
                spectrum[ich] = spectrum[ich] * nbins / specWidth[ich]
            else:
                spectrum[ich] = spectrum[ich] * nbins / specWidth

    #Count rate correction
    rateFactor = 1.
    if not rateStyle == '':
        if isinstance(rateAll, list):
            rateFactor = []
            for ich in range(4):
                rateFactor.append(rateAll[ich] /  float(len(amp[ich])))
        else:
            countAll = 0.0
            for ich in range(4):
                countAll += float(len(amp[ich]))
            #Scale the spectrum by total count rate/total count, so that the integral is total count rate
            rateFactor = rateAll / countAll

    fig = plt.figure(figsize = (12, 8))
    #Single channel plots
    if singlech:
        ich = channel
        gs = gridspec.GridSpec(1, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        ax = fig.add_subplot(gs[0])
        if rateStyle == '':
            #If no count rate correction, the count rate will be calculated as count/total time
            plt.step(x, spectrum / time, where = 'mid', label = 'raw data', zorder = 1)
        else:
            if isinstance(rateAll, list):
                plt.step(x, spectrum * rateFactor[ich], where = 'mid', label = 'raw data', zorder = 1)
            else:
                plt.step(x, spectrum * rateFactor, where = 'mid', label = 'raw data', zorder = 1)
        if parameters.plotSpecLogScale:
            if ecCorr:
                ax.set_xlim([parameters.logScaleBegin, parameters.energyPlotUpper])
            else:
                ax.set_xlim([parameters.logScaleBegin, parameters.adcMax])
            ax.set_xscale('log')
            ax.set_yscale('log')
        else:
            if ecCorr:
                ax.set_xlim([0., parameters.energyPlotUpper])
            else:
                ax.set_xlim([0., parameters.adcMax])
        titleAppend = ''
        if parameters.displayUTC and utc is not None:
            titleAppend = '\n' + convertTimestampToStr(min(utc)) + ' - ' + convertTimestampToStr(max(utc))
        ax.set_title('Spectrum of raw data from ' + filename + titleAppend)
        if ecCorr:
            ax.set_xlabel('energy/keV')
        else:
            ax.set_xlabel('ADC/channel')
        ax.set_ylabel('count rate/cps')
        ax.legend(loc = 0)
        ax.grid()
    #Multiple channel plots
    else:
        gs = gridspec.GridSpec(4, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        for ich in range(4):
            ax = fig.add_subplot(gs[ich])
            if rateStyle == '':
                #If no count rate correction, the count rate will be calculated as count/total time
                plt.step(x[ich], spectrum[ich] / time, where = 'mid', label = 'raw data', zorder = 1)
            else:
                if isinstance(rateAll, list):
                    plt.step(x[ich], spectrum[ich] * rateFactor[ich], where = 'mid', label = 'raw data', zorder = 1)
                else:
                    plt.step(x[ich], spectrum[ich] * rateFactor, where = 'mid', label = 'raw data', zorder = 1)
            if parameters.plotSpecLogScale:
                if ecCorr:
                    ax.set_xlim([parameters.logScaleBegin, parameters.energyPlotUpper])
                else:
                    ax.set_xlim([parameters.logScaleBegin, parameters.adcMax])
                ax.set_xscale('log')
                ax.set_yscale('log')
            else:
                if ecCorr:
                    ax.set_xlim([0., parameters.energyPlotUpper])
                else:
                    ax.set_xlim([0., parameters.adcMax])
            if ich == 0:
                titleAppend = ''
                if parameters.displayUTC and utc is not None:
                    titleAppend = '\n' + convertTimestampToStr(min(utc)) + ' - ' + convertTimestampToStr(max(utc))
                ax.set_title('Spectrum of raw data from ' + filename + titleAppend)
            if ecCorr:
                ax.set_xlabel('energy/keV')
            else:
                ax.set_xlabel('ADC/channel')
            ax.set_ylabel('count rate/cps')
            ax.legend(loc = 0)
            ax.grid()
    if parameters.saveFigNoPlot:
        fig.savefig(saveFigPath + 'spectrum_raw_' + filenameNoAppend + '.png')
        plt.close(fig)
    else:
        plt.show()
    return

def plotLightCurve(filename, uscountEvt, timeStep = 1.0, allCh = True, sciNum = [], curNum = -1, utc = None):

    """
    Function for plotting the time variance curve(light curve) of event counts

    Parameters
    ----------
    filename : str,
        name of the input file\n
    uscountEvt : array of array or array-like,
        event uscount of all 4 channels\n
    timeStep : float, optional
        time step for plotting the variation curve (or bin width of the event time histogram)\n
    allCh : boolean, optional
        indicates plot style, plot light curve of sum of all channels if True, else plot the curve of four channels seperately\n
    sciNum : array of array or array-like, optional
        science data run number of all 4 channels, designed for in-orbit data quick-check\n
    curNum : int, optional
        current scan number for multiple scans data, -1 for single scan\n
    utc : array-like, optional
        UTC readout value, to specify time period of each run for in-orbit data quick-view\n

    Returns
    ----------
    Nothing
    """

    #Save plot option
    filenameNoAppend = filename[0:filename.find(filename.split('.')[-1]) - 1]
    lightSavePath = parameters.lightSaveRootPath + filenameNoAppend + parameters.lightSaveDir
    if parameters.saveFigNoPlot:
        if not os.path.isdir(lightSavePath):
            os.makedirs(lightSavePath)
    #All-channels event plot
    if allCh:
        #Split runs
        sciNumAll = []
        uscountEvtAll = []
        #Merge data from all 4 channels
        for ich in range(4):
            uscountEvtAll += list(uscountEvt[ich])
            if len(sciNum) == 4 and len(sciNum[ich]) > 0:
                sciNumAll += list(sciNum[ich])
        #Split the data according to science runs
        if len(sciNumAll) > 0:
            uscountEvtAll = np.array(uscountEvtAll)
            sciNumAll = np.array(sciNumAll)
            totalRuns = max(sciNumAll) + 1
            #Plot data for each science run seperately
            for irun in range(totalRuns):
                qRun = sciNumAll == irun
                if len(sciNumAll[qRun]) == 0:
                    #Skip empty runs
                    continue
                #If split length specified, double-split the science run data according to the split length
                if parameters.splitRunsLen > 0:
                    curruscountEvtAll = np.array(uscountEvtAll[qRun])
                    splitcount= 0
                    currSplit = 0
                    #Do plot for each split
                    while currSplit < max(curruscountEvtAll):
                        qSplit = (curruscountEvtAll >= currSplit) * (curruscountEvtAll < currSplit + parameters.splitRunsLen)
                        currplot = curruscountEvtAll[qSplit]
                        #UTC values for plot
                        currutc = []
                        if parameters.displayUTC and utc is not None:
                            currutc = utc[qSplit]
                        if len(currplot) > 0:
                            #Only plot non-empty datasets
                            specusc, xusc = np.histogram(currplot, bins = int((np.max(currplot) - np.min(currplot) + 10) / timeStep), range = (np.min(currplot) - 5, \
                                np.max(currplot) + 5))
                            xusc = (xusc[:-1] + xusc[1:]) / 2
                            specusc = specusc / (((np.max(currplot) - np.min(currplot) + 10)) / float(int((np.max(currplot) - np.min(currplot) + 10) / timeStep)) * \
                                np.ones(len(specusc))) #Divide bin count by bin width
                            fig = plt.figure(figsize = (12, 8))
                            plt.step(xusc, specusc)
                            plt.xlim([np.min(currplot) - 5, np.max(currplot) + 5])
                            plt.xlabel('time/s')
                            plt.ylabel('count rate/cps')
                            titleAppend = ''
                            if parameters.displayUTC and len(currutc) > 0:
                                titleAppend = '\n' + convertTimestampToStr(min(currutc)) + ' - ' + convertTimestampToStr(max(currutc))
                            plt.title('Count rate variation of file ' + filename + ', run #' + str(irun) + '-' + str(splitcount) + titleAppend)
                            plt.grid()
                            if parameters.saveFigNoPlot:
                                fig.savefig(lightSavePath + 'light_curve_' + str(irun) + '_' + str(splitcount) + '_' + filenameNoAppend + '.png')
                                plt.close(fig)
                            else:
                                plt.show()
                        currSplit += parameters.splitRunsLen
                        splitcount += 1
                #If split length not specified, plot data of each science run
                else:
                    specusc, xusc = np.histogram(uscountEvtAll[qRun], bins = int(1.2 * np.max(uscountEvtAll[qRun]) / timeStep), range = (- 0.1 * \
                        np.max(uscountEvtAll[qRun]), 1.1 * np.max(uscountEvtAll[qRun])))
                    xusc = (xusc[:-1] + xusc[1:]) / 2
                    specusc = specusc / (1.2 * np.max(uscountEvtAll[qRun]) / float(int(1.2 * np.max(uscountEvtAll[qRun]) / timeStep)) * \
                        np.ones(len(specusc))) #Divide bin count by bin width
                    fig = plt.figure(figsize = (12, 8))
                    plt.step(xusc, specusc)
                    plt.xlabel('time/s')
                    plt.ylabel('count rate/cps')
                    titleAppend = ''
                    if parameters.displayUTC and utc is not None:
                        titleAppend = '\n' + convertTimestampToStr(min(utc)) + ' - ' + convertTimestampToStr(max(utc))
                    plt.title('Count rate variation of file ' + filename + ', run #' + str(irun) + titleAppend)
                    plt.grid()
                    if parameters.saveFigNoPlot:
                        fig.savefig(lightSavePath + 'light_curve_' + str(irun) + '_' + filenameNoAppend + '.png')
                        plt.close(fig)
                    else:
                        plt.show()
        #If not split, plot all data
        else:
            specusc, xusc = np.histogram(uscountEvtAll, bins = int(1.2 * np.max(uscountEvtAll) / timeStep), range = (- 0.1 * np.max(uscountEvtAll), \
                1.1 * np.max(uscountEvtAll)))
            xusc = (xusc[:-1] + xusc[1:]) / 2
            specusc = specusc / (1.2 * np.max(uscountEvtAll) / float(int(1.2 * np.max(uscountEvtAll) / timeStep)) * \
                np.ones(len(specusc))) #Divide bin count by bin width
            fig = plt.figure(figsize = (12, 8))
            plt.step(xusc, specusc)
            plt.xlabel('time/s')
            plt.ylabel('count rate/cps')
            titleAppend = ''
            if parameters.displayUTC and utc is not None:
                titleAppend = '\n' + convertTimestampToStr(min(utc)) + ' - ' + convertTimestampToStr(max(utc))
            if curNum < 0:
                plt.title('Count rate variation of file ' + filename + titleAppend)
            else:
                plt.title('Count rate variation of file ' + filename + ', scan #' + str(curNum) + titleAppend)
            plt.grid()
            if parameters.saveFigNoPlot:
                if curNum < 0:
                    fig.savefig(lightSavePath + 'light_curve_' + filenameNoAppend + '.png')
                else:
                    fig.savefig(lightSavePath + 'light_curve_' + str(curNum) + '_' + filenameNoAppend + '.png')
                plt.close(fig)
            else:
                plt.show()
    #Plot events of different channels seperately
    else:
        #Split runs
        if len(sciNum) == 4:
            sciNumAll = []
            #Merge science run numbers to get total run numbers
            for ich in range(4):
                if len(sciNum) == 4 and len(sciNum[ich]) > 0:
                    sciNumAll += list(sciNum[ich])
            #Split the data according to science runs
            if len(sciNumAll) > 0:
                sciNumAll = np.array(sciNumAll)
                totalRuns = max(sciNumAll) + 1
                #Plot data for each science run seperately
                for irun in range(totalRuns):
                    qRun = []
                    for ich in range(4):
                        qRun.append(np.array(sciNum[ich]) == irun)
                    #If split length specified, double-split the science run data according to the split length
                    if parameters.splitRunsLen > 0:
                        curruscountEvt = []
                        runMax = -1
                        #Get maximum uscount value in current run
                        for ich in range(4):
                            curruscountEvt.append(np.array(uscountEvt[ich][qRun[ich]]))
                            if len(curruscountEvt[ich]) > 0 and max(curruscountEvt[ich]) > runMax:
                                runMax = max(curruscountEvt[ich])
                        curruscountEvt = np.array(curruscountEvt)
                        splitcount= 0
                        currSplit = 0
                        #Do plot for each split
                        while currSplit < runMax:
                            fig = plt.figure(figsize = (12, 8))
                            gs = gridspec.GridSpec(4, 1, wspace=0.5, hspace=0.2, left=0.13, right=0.95)
                            iFirst = -1
                            currplot = []
                            #Find first non-empty channel
                            for ich in range(4):
                                if len(curruscountEvt[ich]) > 0:
                                    qSplit = (np.array(curruscountEvt[ich]) >= currSplit) * (np.array(curruscountEvt[ich]) < currSplit + parameters.splitRunsLen)
                                    currplot.append(np.array(curruscountEvt[ich])[qSplit])
                                else:
                                    #Skip empty channels
                                    currplot.append([])
                            for ich in range(4):
                                if len(currplot[ich]) > 0:
                                    #Find first non-empty channel
                                    iFirst = ich
                                    break
                            #UTC values for plot
                            currutc = []
                            if parameters.displayUTC and utc is not None:
                                currutc = utc[qSplit]
                            #Do plot if not all channels are empty
                            if iFirst > -1:
                                tmin = min(currplot[iFirst])
                                tmax = max(currplot[iFirst])
                                #Get plot range
                                for ich in range(4):
                                    if len(currplot[ich]) == 0:
                                        #Skip empty channels
                                        continue
                                    if max(currplot[ich]) > tmax:
                                        tmax = max(currplot[ich])
                                    if min(currplot[ich]) < tmin:
                                        tmin = min(currplot[ich])
                                ax = None
                                for ich in range(4):
                                    if len(np.array(currplot[ich])) == 0:
                                        #Only plot non-empty datasets
                                        continue
                                    ax = fig.add_subplot(gs[ich])
                                    specusc, xusc = np.histogram(currplot[ich], bins = int((np.max(currplot[ich]) - np.min(currplot[ich]) + 10) / timeStep), range = \
                                        (np.min(currplot[ich]) - 5, np.max(currplot[ich]) + 5))
                                    xusc = (xusc[:-1] + xusc[1:]) / 2
                                    specusc = specusc / ((np.max(currplot[ich]) - np.min(currplot[ich]) + 10) / float(int((np.max(currplot[ich]) - np.min(currplot[ich]) + 10) \
                                        / timeStep)) * np.ones(len(specusc))) #Divide bin count by bin width
                                    ax.step(xusc, specusc)
                                    ax.set_ylabel('count rate/cps')
                                    if ich == 0:
                                        titleAppend = ''
                                        if parameters.displayUTC and len(currutc) > 0:
                                            titleAppend = '\n' + convertTimestampToStr(min(currutc)) + ' - ' + convertTimestampToStr(max(currutc))
                                        ax.set_title('Count rate variation of file ' + filename + ', run#' + str(irun) + '-' + str(splitcount) + titleAppend)
                                    ax.grid()
                                    ax.set_xlim([tmin - 5, tmax + 5])
                                try:
                                    ax.set_xlabel('time/s')
                                    if parameters.saveFigNoPlot:
                                        fig.savefig(lightSavePath + 'light_curve_' + str(irun) + '_' + str(splitcount) + '_' + filenameNoAppend + '.png')
                                    else:
                                        plt.show()
                                except:
                                    #Skip empty channels
                                    pass
                            if parameters.saveFigNoPlot:
                                plt.close(fig)
                            currSplit += parameters.splitRunsLen
                            splitcount += 1
                    #If split length not specified, plot data of each science run
                    else:
                        fig = plt.figure(figsize = (12, 8))
                        gs = gridspec.GridSpec(4, 1, wspace=0.5, hspace=0.2, left=0.13, right=0.95)
                        iFirst = -1
                        #Find first non-empty channel
                        for ich in range(4):
                            if len(np.array(uscountEvt[ich])[qRun[ich]]) > 0:
                                iFirst = ich
                                break
                        #Do plot if not all channels are empty
                        if iFirst > -1:
                            tmin = min(np.array(uscountEvt[iFirst])[qRun[iFirst]])
                            tmax = max(np.array(uscountEvt[iFirst])[qRun[iFirst]])
                            for ich in range(4):
                                if len(np.array(uscountEvt[ich])[qRun[ich]]) == 0:
                                    #Skip empty channels
                                    continue
                                if max(np.array(uscountEvt[ich])[qRun[ich]]) > tmax:
                                    tmax = max(np.array(uscountEvt[ich])[qRun[ich]])
                                if min(np.array(uscountEvt[ich])[qRun[ich]]) < tmin:
                                    tmin = min(np.array(uscountEvt[ich])[qRun[ich]])
                            ax = None
                            for ich in range(4):
                                if len(np.array(uscountEvt[ich])[qRun[ich]]) == 0:
                                    #Only plot non-empty datasets
                                    continue
                                ax = fig.add_subplot(gs[ich])
                                specusc, xusc = np.histogram(np.array(uscountEvt[ich])[qRun[ich]], bins = int(1.2 * np.max(np.array(uscountEvt[ich])[qRun[ich]]) \
                                    / timeStep), range = (- 0.1 * np.max(np.array(uscountEvt[ich])[qRun[ich]]), 1.1 * np.max(np.array(uscountEvt[ich])[qRun[ich]])))
                                xusc = (xusc[:-1] + xusc[1:]) / 2
                                specusc = specusc / (1.2 * np.max(np.array(uscountEvt[ich])[qRun[ich]]) / float(int(1.2 * np.max(np.array(uscountEvt[ich])[qRun[ich]]) \
                                    / timeStep)) * np.ones(len(specusc))) #Divide bin count by bin width
                                ax.step(xusc, specusc)
                                ax.set_ylabel('count rate/cps')
                                if ich == 0:
                                    titleAppend = ''
                                    if parameters.displayUTC and utc is not None:
                                        titleAppend = '\n' + convertTimestampToStr(min(utc)) + ' - ' + convertTimestampToStr(max(utc))
                                    ax.set_title('Count rate variation of file ' + filename + ', run#' + str(irun) + titleAppend)
                                ax.grid()
                                ax.set_xlim([tmin, tmax])
                            try:
                                ax.set_xlabel('time/s')
                                if parameters.saveFigNoPlot:
                                    fig.savefig(lightSavePath + 'light_curve_' + str(irun) + '_' + filenameNoAppend + '.png')
                                else:
                                    plt.show()
                            except:
                                #Skip empty channels
                                pass
                        if parameters.saveFigNoPlot:
                            plt.close(fig)
                return
        
        #Plot all data, no science run split
        fig = plt.figure(figsize = (12, 8))
        gs = gridspec.GridSpec(4, 1, wspace=0.5, hspace=0.2, left=0.13, right=0.95)
        tmin = min(uscountEvt[0])
        tmax = max(uscountEvt[0])
        #Get plot range
        for ich in range(1, 4):
            if max(uscountEvt[ich]) > tmax:
                tmax = max(uscountEvt[ich])
            if min(uscountEvt[ich]) < tmin:
                tmin = min(uscountEvt[ich])
        ax = None
        for ich in range(4):
            ax = fig.add_subplot(gs[ich])
            specusc, xusc = np.histogram(uscountEvt[ich], bins = int(1.2 * np.max(uscountEvt[ich]) / timeStep), range = (- 0.1 * \
                np.max(uscountEvt[ich]), 1.1 * np.max(uscountEvt[ich])))
            xusc = (xusc[:-1] + xusc[1:]) / 2
            specusc = specusc / (1.2 * np.max(uscountEvt[ich]) / float(int(1.2 * np.max(uscountEvt[ich]) / timeStep)) * \
                np.ones(len(specusc))) #Divide bin count by bin width
            ax.step(xusc, specusc)
            ax.set_ylabel('count rate/cps')
            if ich == 0:
                titleAppend = ''
                if parameters.displayUTC and utc is not None:
                    titleAppend = '\n' + convertTimestampToStr(min(utc)) + ' - ' + convertTimestampToStr(max(utc))
                if curNum < 0:
                    ax.set_title('Count rate variation of file ' + filename + titleAppend)
                else:
                    ax.set_title('Count rate variation of file ' + filename + ', scan #' + str(curNum) + titleAppend)
            ax.grid()
            ax.set_xlim([tmin, tmax])
        ax.set_xlabel('time/s')
        if parameters.saveFigNoPlot:
            if curNum < 0:
                fig.savefig(lightSavePath + 'light_curve_' + filenameNoAppend + '.png')
            else:
                fig.savefig(lightSavePath + 'light_curve_' + str(curNum) + '_' + filenameNoAppend + '.png')
            plt.close(fig)
        else:
            plt.show()
    return

#*************************************************************************************************************************************************************
#************************************************Basic calibration (temperature, bias and EC) part************************************************
#*************************************************************************************************************************************************************
#TBD: add a input parameter to read the bin factor/fit range from input files
def getBiasnbinsFactor(scanNum, source, scanRange = []):
    
    """
    NOTE: Used for old Na-22 bias scan(multiple scan) data only
    ----------

    Function for calculating the binning factor for bias response process

    Parameters
    ----------
    scanNum : int,
        scan number\n
    source : str, 
        name of the source, currently supporting only `Am241` and `Na22`\n
    scanRange : list, optional
        containing the accepted scan range for multiple scans\n

    Returns
    ----------
    int,
        corresponding binning factor(in nbins/512)\n
    """

    if source == 'Na22':
        biasnbinsFactor = [8, 4, 2, 1]
        biasnbinsRange = [1, 5, 12, 18]
    else:
        biasnbinsFactor = [2, 1]
        biasnbinsRange = [10, 17]

    num = scanNum
    if not len(scanRange) == 0:
        num += scanRange[0]
    ind = 0
    while ind < len(biasnbinsRange) - 1:
        if num >= biasnbinsRange[ind] and num < biasnbinsRange[ind + 1]:
            break
        ind += 1
    return biasnbinsFactor[ind]

def getBiasFitRange(scanNum, source, corr = False):
    
    """
    NOTE: Used for old Na-22 bias scan(multiple scan) data only
    ----------

    Function for calculating the fit range for bias response process

    Parameters
    ----------
    scanNum : int,
        scan number\n
    source : str, 
        name of the source, currently supporting only `Am241` and `Na22`\n
    corr : boolean, optional
        indicates whether the data is temperature-calibrated\n

    Returns
    ----------
    list,
        corresponding fit range\n
    """

    if corr:
        biasRange = [1, 3, 7, 14]
        biasFitRange = [[[1000, 3250], [1250, 3800], [1000, 3250], [1200, 3500]], [[2500, 7000], [3000, 7500], [2500, 7000], [2500, 7000]], \
            [[5000, 15000], [6000, 16000], [5000, 15000], [5500, 15500]], [[12000, 26000], [13000, 29000], [12000, 28000], [13000, 29000]]]
        num = scanNum + 1
        ind = 0
        while ind < len(biasFitRange) - 1:
            if num >= biasRange[ind] and num < biasRange[ind + 1]:
                break
            ind += 1
        return biasFitRange[ind]
    else:
        biasFitRange = []
        if source == 'Am241':
            biasFitRange = [[], [], [], [], [], [], [[600, 800], [640, 900], [580, 800], [630, 900]], [[600, 800], [650, 1000], [580, 900], [630, 900]], \
                [[600, 1000], [650, 1100], [600, 1000], [650, 1100]], [[600, 1100], [650, 1300], [600, 1200], [650, 1200]], \
                [[600, 1300], [650, 1400], [600, 1300], [650, 1400]], [[600, 1500], [750, 1750], [600, 1600], [650, 1600]], \
                [[700, 1700], [800, 1800], [700, 1700], [700, 1700]], [], [], [], [[1000, 2300], [1300, 2600], [1000, 2500], [1000, 2600]], \
                [[1100, 2600], [1400, 3000], [1100, 2600], [1200, 2800]], [], [], []]
        elif source == 'Na22':
            biasFitRange = [[[1000, 1750], [1250, 2250], [1000, 1750], [1200, 2000]],
                [[1400, 2400], [1750, 3000], [1500, 2500], [1600, 2600]],
                [[2000, 3250], [2600, 3800], [2000, 3250], [2250, 3500]],
                [[2500, 4000], [3000, 5000], [2500, 4000], [2500, 4500]],
                [[3000, 5000], [3500, 5500], [3000, 5000], [3000, 5000]],
                [[3500, 6000], [4000, 6500], [3500, 6000], [4000, 6500]],
                [[4500, 7000], [5000, 7500], [4500, 7000], [4500, 7000]],
                [[5000, 8000], [6000, 8500], [5000, 8000], [5500, 8500]],
                [[6000, 9000], [6500, 10000], [6000, 9000], [6500, 9500]],
                [[7000, 10500], [8000, 12000], [7000, 12000], [7500, 12000]],
                [[8000, 12000], [8500, 13000], [8000, 13000], [8500, 13000]],
                [[9000, 14000], [11000, 15000], [9000, 14500], [9500, 15000]],
                [[10000, 15000], [11000, 16000], [10000, 16000], [11500, 16000]],
                [[11000, 17000], [12000, 18000], [11000, 18000], [11500, 18000]],
                [[12000, 18000], [13000, 20000], [13000, 20000], [13000, 20000]],
                [[13000, 20000], [14500, 21000], [14000, 21000], [14500, 22000]],
                [[14000, 22000], [15000, 23000], [14000, 23000], [14000, 24000]],
                [[16000, 24000], [17000, 25000], [18000, 25000], [17000, 26000]],
                [[17000, 26000], [20000, 27000], [17000, 27000], [18000, 28000]],
                [[19000, 28000], [20000, 29000], [18000, 29000], [19000, 30000]],
                [[20000, 30000], [24000, 32000], [23000, 32000], [24000, 33000]]]
        return biasFitRange[scanNum]

def ecCorrection(amp, singlech = False, channel = -1, ecFiles = parameters.ecConfig, dataCorr = True):
    #TBD: Add error to corrected energy values
    """
    NOTE: To be implemented in v1.0.0 complete
    ----------

    Function for calculating the EC calibrated values of the spectrum's x-axis(energy)

    Parameters
    ----------
    amp : array of arrays or array-like,
        the input ADC values\n
    singlech : boolean, optional
        indicates whether the input data is single-channeled\n
    channel : int, optional
        the channel number in range [0-3]\n
    ecFiles : list, optional
        list containing names of files with EC correction coeficients\n
    dataCorr : boolean, optional
        indicates whether the input is amplitude data, False when input is spectrum fit range to be corrected

    Returns
    ----------
    corrAmp : array of arrays or array-like,
        the calibrated ADC values(converted from ADC channel to energy)\n
    
    Raises
    ----------
    Exception
        when error reading EC coefficient file, or channel number of EC coefficient file does not match current channel
    """

    corrAmp = copy(amp)
    if singlech:
        #Read EC coefficients from config file
        try:
            with open(ecFiles[channel], 'r') as fin:
                ecCoef = json.load(fin)
                if ecCoef['channel'] != channel:
                    raise Exception()
        except:
            raise Exception('tempBiasCorrection: error reading config file \"' + ecFiles[channel] + '\"')
        #EC correction for amplitude data
        if dataCorr:
            #Low energy region
            qLow = np.array(amp) <= ecCoef['EC_bound_low']
            corrAmp[qLow] = quadFunction(ecCoef['EC_low'], np.array(amp)[qLow])
            #Medium energy region
            qMid = (np.array(amp) > ecCoef['EC_bound_low']) * (np.array(amp) < ecCoef['EC_bound_high'])
            corrAmp[qMid] = linearFunction(ecCoef['EC_mid'], np.array(amp)[qMid])
            #High energy region
            qHigh = np.array(amp) >= ecCoef['EC_bound_high']
            corrAmp[qHigh] = quadFunction(ecCoef['EC_high'], np.array(amp)[qHigh])
        #EC correction for spectrum fit range
        else:
            if len(corrAmp) == 0:
                pass
            if amp[0] <= ecCoef['EC_bound_low']:
                corrAmp[0] = quadFunction(ecCoef['EC_low'], amp[0])
            elif amp[0] <= ecCoef['EC_bound_high']:
                corrAmp[0] = linearFunction(ecCoef['EC_mid'], amp[0])
            else:
                corrAmp[0] = quadFunction(ecCoef['EC_high'], amp[0])
            if amp[1] <= ecCoef['EC_bound_low']:
                corrAmp[1] = quadFunction(ecCoef['EC_low'], amp[1])
            elif amp[1] <= ecCoef['EC_bound_high']:
                corrAmp[1] = linearFunction(ecCoef['EC_mid'], amp[1])
            else:
                corrAmp[1] = quadFunction(ecCoef['EC_high'], amp[1])
    else:
        for ich in range(4):
            #Read EC coefficients from config file
            try:
                with open(ecFiles[ich], 'r') as fin:
                    ecCoef = json.load(fin)
                    if ecCoef['channel'] != ich:
                        raise Exception()
            except:
                raise Exception('tempBiasCorrection: error reading config file \"' + ecFiles[ich] + '\"')
            #EC correction for amplitude data
            if dataCorr:
                corrAmp[ich] = np.array(corrAmp[ich])
                #Low energy region
                qLow = np.array(amp[ich]) <= ecCoef['EC_bound_low']
                corrAmp[ich][qLow] = quadFunction(ecCoef['EC_low'], np.array(amp[ich])[qLow])
                #Medium energy region
                qMid = (np.array(amp[ich]) > ecCoef['EC_bound_low']) * (np.array(amp[ich]) < ecCoef['EC_bound_high'])
                corrAmp[ich][qMid] = linearFunction(ecCoef['EC_mid'], np.array(amp[ich])[qMid])
                #High energy region
                qHigh = np.array(amp[ich]) >= ecCoef['EC_bound_high']
                corrAmp[ich][qHigh] = quadFunction(ecCoef['EC_high'], np.array(amp[ich])[qHigh])
            #EC correction for spectrum fit range
            else:
                if len(corrAmp[ich]) == 0:
                    pass
                if amp[ich][0] <= ecCoef['EC_bound_low']:
                    corrAmp[ich][0] = quadFunction(ecCoef['EC_low'], amp[ich][0])
                elif amp[ich][0] <= ecCoef['EC_bound_high']:
                    corrAmp[ich][0] = linearFunction(ecCoef['EC_mid'], amp[ich][0])
                else:
                    corrAmp[ich][0] = quadFunction(ecCoef['EC_high'], amp[ich][0])
                if amp[ich][1] <= ecCoef['EC_bound_low']:
                    corrAmp[ich][1] = quadFunction(ecCoef['EC_low'], amp[ich][1])
                elif amp[ich][1] <= ecCoef['EC_bound_high']:
                    corrAmp[ich][1] = linearFunction(ecCoef['EC_mid'], amp[ich][1])
                else:
                    corrAmp[ich][1] = quadFunction(ecCoef['EC_high'], amp[ich][1])

    return corrAmp

def tempBiasCorrection(temp, bias, doCorr = True, coefFile = parameters.tempBiasConfig):

    """
    Function for calculating correction factor of temperature and bias

    Parameters
    ----------
    temp : array of arrays,
        temperature of SiPM\n
    bias : array of arrays,
        SiPM bias\n
    doCorr : boolean, optional
        indicates whether the temperature-bias correction will be done, to avoid warning info output\n
    coefFile : str, optional
        the file for correction coeficients, DO NOT CHANGE UNLESS NECESSARY\n

    Returns
    ----------
    corrFactor : list
        temp-bias correction factors of all 4 channels, in the form of [float, float, float, float]
    corrErr : list
        error of correction factors, in the same form as corrFactor
    
    Raises
    ----------
    Exception
        when error reading temp-bias coefficient file
    """

    def quad2dFunctionInternal(param, x, y, xdiff = 0, ydiff = 0):

        """
        Internal function to calculate value of 2-d quadratic function and its (1st) derivative,
        with the formula being\n
        `f(x, y) = p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ^ 2 + p02 * y ^ 2 + p21 * x ^ 2 * y + p12 * x * y ^ 2 + p22 * x ^ 2 * y ^ 2 + p30 + x ^ 3 + \
            p31 * x * 3 * y + p40 * x ^ 4`

        Parameters
        ----------
        param : list,
            parameters of quadratic function, in the form of [p00, p10, p01, p11, p20, p02, p21, p12, p22, p30, p31, p40]\n
        x : float or array-like,
            value(s) of input x\n
        y : float or array-like,
            value(s) of input y\n
        xdiff : int, optional
            partial derivative order of x, default 0 for `f(x, y)`, and 1 for `df/dx`\n
        ydiff : int, optional
            partial derivative order of y, default 0 for `f(x, y)`, and 1 for `df/dy`\n

        Returns
        ----------
        quad : float or array-like,
            values of quadratic function corresponding to input x\n

        Raises
        ----------
        Exception
            when xdiff or ydiff is not 0 or 1
        """

        quad = None
        if xdiff == 0:
            if ydiff == 0:
                quad = param[0] + param[1] * x + param[2] * y + param[3] * x * y + param[4] * x ** 2 + param[5] * y ** 2 + param[6] * x ** 2 * y + \
                    param[7] * x * y ** 2 + param[8] * x ** 2 * y ** 2 + param[9] * x ** 3 + param[10] * x ** 3 * y + param[11] * x ** 4
            elif ydiff == 1:
                quad = param[2] + param[3] * x + 2 * param[5] * y + param[6] * x ** 2 + 2 * param[7] * x * y + 2 * param[8] * x ** 2 * y  + param[10] * x ** 3
            else:
                raise Exception('quadFunctionInternal: illegal derivative order ' + str(ydiff))
        elif xdiff == 1:
            if ydiff == 0:
                quad = param[1] + param[3] * y + 2 * param[4] * x + 2 * param[6] * x * y + param[7] * y ** 2 + 2 * param[8] * x * y ** 2 + 3 * param[9] * x ** 2 + \
                    3 * param[10] * x ** 2 * y + 4 * param[11] * x ** 3
            elif ydiff == 1:
                quad = param[3] + 2 * param[6] * x + 2 * param[7] * y + 4 * param[8] * x * y + 3 * param[10] * x ** 2
            else:
                raise Exception('quadFunctionInternal: illegal derivative order ' + str(ydiff))
        else:
            raise Exception('quadFunctionInternal: illegal derivative order ' + str(xdiff))
        return quad

    if not doCorr:
        return [1.0, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0]

    #Standard temperature and bias
    tempStandard = parameters.tempStandard
    biasStandard = parameters.biasStandard
    corrFactor = []
    corrErr = []

    for ich in range(4):
        tempAvg = np.average(temp[ich])
        tempStd = np.sqrt(np.std(temp[ich]) ** 2 + parameters.tempErrSys ** 2)
        biasAvg = np.average(bias[ich])
        biasStd = np.std(bias[ich])
        #2-D(correlated) temp-bias correction
        #Read config file
        try:
            with open(coefFile, 'r') as fin:
                tempBiasCoef = json.load(fin)
        except:
            raise Exception('tempBiasCorrection: error reading config file \"' + coefFile + '\"')

        #Coefficients of 2-D(correlated) temp-bias correction
        #2-D temp-bias response: p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ** 2 + p 02 * y ** 2 + p21 * x ** 2 * y +
        #p12 * x * y ** 2 + p22 * x ** 2 * y ** 2 + p30 * x ** 3 + p31 * x ** 3 * y + p 40 ** 4
        p00 = [tempBiasCoef[ich]['p00'] for ich in range(4)]
        p00Err = [tempBiasCoef[ich]['p00_err'] for ich in range(4)]
        p10 = [tempBiasCoef[ich]['p10'] for ich in range(4)]
        p10Err = [tempBiasCoef[ich]['p10_err'] for ich in range(4)]
        p01 = [tempBiasCoef[ich]['p01'] for ich in range(4)]
        p01Err = [tempBiasCoef[ich]['p01_err'] for ich in range(4)]
        p11 = [tempBiasCoef[ich]['p11'] for ich in range(4)]
        p11Err = [tempBiasCoef[ich]['p11_err'] for ich in range(4)]
        p20 = [tempBiasCoef[ich]['p20'] for ich in range(4)]
        p20Err = [tempBiasCoef[ich]['p20_err'] for ich in range(4)]
        p02 = [tempBiasCoef[ich]['p02'] for ich in range(4)]
        p02Err = [tempBiasCoef[ich]['p02_err'] for ich in range(4)]
        p21 = [tempBiasCoef[ich]['p21'] for ich in range(4)]
        p21Err = [tempBiasCoef[ich]['p21_err'] for ich in range(4)]
        p12 = [tempBiasCoef[ich]['p12'] for ich in range(4)]
        p12Err = [tempBiasCoef[ich]['p12_err'] for ich in range(4)]
        p22 = [tempBiasCoef[ich]['p22'] for ich in range(4)]
        p22Err = [tempBiasCoef[ich]['p22_err'] for ich in range(4)]
        p30 = [tempBiasCoef[ich]['p30'] for ich in range(4)]
        p30Err = [tempBiasCoef[ich]['p30_err'] for ich in range(4)]
        p31 = [tempBiasCoef[ich]['p31'] for ich in range(4)]
        p31Err = [tempBiasCoef[ich]['p31_err'] for ich in range(4)]
        p40 = [tempBiasCoef[ich]['p40'] for ich in range(4)]
        p40Err = [tempBiasCoef[ich]['p40_err'] for ich in range(4)]

        param = [p00[ich], p10[ich], p01[ich], p11[ich], p20[ich], p02[ich], p21[ich], p12[ich], p22[ich], p30[ich], p31[ich], p40[ich]]
        paramErr = [p00Err[ich], p10Err[ich], p01Err[ich], p11Err[ich], p20Err[ich], p02Err[ich], p21Err[ich], p12Err[ich], p22Err[ich], p30Err[ich], \
            p31Err[ich], p40Err[ich]]
        corrFactor.append(quad2dFunctionInternal(param, tempStandard, biasStandard) / quad2dFunctionInternal(param, tempAvg, biasAvg))
        tempStdTerm = quad2dFunctionInternal(param, tempStandard, biasStandard) * quad2dFunctionInternal(param, tempAvg, biasAvg, \
            xdiff = 1) /quad2dFunctionInternal(param, tempAvg, biasAvg) ** 2
        biasStdTerm = quad2dFunctionInternal(param, tempStandard, biasStandard) * quad2dFunctionInternal(param, tempAvg, biasAvg, \
            ydiff = 1) / quad2dFunctionInternal(param, tempAvg, biasAvg) ** 2
        pTermStandard = [1., tempStandard, biasStandard, tempStandard * biasStandard, tempStandard ** 2, biasStandard ** 2, \
            tempStandard ** 2 * biasStandard, tempStandard * biasStandard ** 2, tempStandard ** 2 * biasStandard ** 2, tempStandard ** 3, \
                tempStandard ** 3 * biasStandard, tempStandard ** 4]
        pTermAvg = [1., tempAvg, biasAvg, tempAvg * biasAvg, tempAvg ** 2, biasAvg ** 2, tempAvg ** 2 * biasAvg, tempAvg * biasAvg ** 2, \
            tempAvg ** 2 * biasAvg ** 2, tempAvg ** 3, tempAvg ** 3 * biasAvg, tempAvg ** 4]
        corrErrSum = (tempStdTerm * tempStd) ** 2 + (biasStdTerm * biasStd) ** 2
        for it in range(len(param)):
            corrErrSum += ((pTermStandard[it] * quad2dFunctionInternal(param, tempAvg, biasAvg) - pTermAvg[it] * quad2dFunctionInternal(\
                param, tempStandard, biasStandard)) / quad2dFunctionInternal(param, tempAvg, biasAvg) ** 2 * paramErr[it]) ** 2
        corrErr.append(np.sqrt(corrErrSum))
    
    return corrFactor, corrErr

#*************************************************************************************************************************************************************
#************************************************************Basic fit functions and fit part***********************************************************
#*************************************************************************************************************************************************************

def gehrelsErr(ydata):

    """
    Function for calculating Poissons error of spectrum count (y) data\n
    in the form of: \n
        sqrt(ydata) if ydata >= 5.0
        else 1 + sqrt(ydata + 0.75)

    Parameters
    ----------
    ydata : array-like or float,
        spectrum data (or other Poissions variables)\n

    Returns
    ----------
    yerr : array-like or float,
        corresponding Poissons error\n
    """

    if np.size(ydata) == 1:
        #Float input
        yerr = np.sqrt(ydata) if ydata >= 5.0 else 1.0 + np.sqrt(ydata + 0.75)
    elif np.size(ydata) > 1:
        #Array input
        yerr = np.sqrt(ydata)
        q = np.where(ydata < 5.0)
        yerr[q] = 1.0 + np.sqrt(ydata[q] + 0.75)
    else:
        #Empty input
        yerr = []
    return yerr

#------------------------------------------------------------Spectrum fit functions------------------------------------------------------------

def gaussianFunction(param, x):

    """
    Auxiliary function to calculate the value of the (normalized) gaussian function, used for ODR fit\n
    in the form of\n
    `g(x) = amplitude * exp(-(x - center) ^ 2 / (2 * peak_sigma ^ 2)) / sqrt(2 * pi * sigma ^ 2)`\n

    Parameters
    ----------
    param : list,
        parameters of the gaussian function, in the form of [amplitude, center, sigma]\n
    x : float or array-like,
        input x value(s)\n

    Returns
    ----------
    float or array-like,
        corresponding value(s) of the gaussian function\n
    """

    return param[0] * np.exp(-(x - param[1]) ** 2 / (2 * (param[2] ** 2))) / (param[2] * np.sqrt(2 * np.pi))

def doFitGaussian(xdata, ydata, odr = False, xerror = [], yerror = []):

    """
    Function for fitting single gaussian peak, with no background

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_amplitude': area of the gaussian peak,
            'fit_center': center of the gaussian peak,
            'fit_sigma': standard deviation of the gaussian peak,
            'fit_amplitude_err': error of area of the gaussian peak,
            'fit_center_err': error of center of the gaussian peak,
            'fit_sigma_err': error of standard deviation of the gaussian peak,
        }\n
    with fit function as\n
    `y = fit_amplitude * exp(-(x - fit_center) ^ 2 / (2 * fit_sigma ^ 2)) / sqrt(2 * pi * fit_sigma ^ 2)`\n
    """

    gModel = lmfit.models.GaussianModel(prefix = 'fit_')
    #Initial guess of parameters
    param = gModel.guess(ydata, x = xdata)
    if odr:
        #ODR(Orthogonal distance regression) fit
        if len(xerror) == 0:
            if len(yerror) == 0:
                data = RealData(xdata, ydata)
            else:
                data = RealData(xdata, ydata, sy = yerror)
        else:
            if len(yerror) == 0:
                data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
            else:
                data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
        model = Model(gaussianFunction)
        odrFit = ODR(data, model, [param.valuesdict()['fit_amplitude'], param.valuesdict()['fit_center'], param.valuesdict()['fit_sigma']])
        odrFit.set_job(fit_type = 0)
        result = odrFit.run()
        fitResult = {
                        'fit_amplitude':           result.beta[0],
                        'fit_center':                result.beta[1],
                        'fit_sigma':                 result.beta[2],
                        'fit_amplitude_err':           result.sd_beta[0],
                        'fit_center_err':                result.sd_beta[1],
                        'fit_sigma_err':                 result.sd_beta[2],
            }
    else:
        #Least-squares fit
        if len(yerror) == 0:
            result = gModel.fit(ydata, param, x = xdata)
        else:
            result = gModel.fit(ydata, param, x = xdata, weights = 1. / yerror)
        fitResult = {
                        'fit_amplitude':           result.best_values['fit_amplitude'],
                        'fit_center':                result.best_values['fit_center'],
                        'fit_sigma':                 result.best_values['fit_sigma'],
                        'fit_amplitude_err':           result.params['fit_amplitude'].stderr,
                        'fit_center_err':                result.params['fit_center'].stderr,
                        'fit_sigma_err':                 result.params['fit_sigma'].stderr,
            }
    return fitResult

def peakFunctionLin(param, x):

    """
    Auxiliary function to calculate the value of the peak function with linear background, used for odr fitting\n
    in the form of\n
    `g(x) = amplitude * exp(-(x - center) ^ 2 / (2 * peak_sigma ^ 2)) / sqrt(2 * pi * sigma ^ 2) + a * x + b`\n

    param : list,
        parameters of the gaussian and linear function, in the form of [amplitude, center, sigma, a, b]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value of the peak function, as an gaussian + linear\n
    """

    return param[0] * np.exp(-(x - param[1]) ** 2 / (2 * (param[2] ** 2))) / (param[2] * np.sqrt(2 * np.pi)) + param[3] * x + param[4]

def peakFunctionQuad(param, x):

    """
    Auxiliary function to calculate the value of the peak function with quadratic background, used for ODR fit\n
    in the form of\n
    `g(x) = amplitude * exp(-(x - center) ** 2 / (2 * sigma ** 2)) / sqrt(2 * pi * sigma ^ 2) + a * x ** 2 + b * x + c`\n

    param : list,
        parameters of the gaussian and quadratic function, in the form of [amplitude, center, sigma, a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value of the peak function, as an gaussian + quadratic\n
    """

    return param[0] * np.exp(-(x - param[1]) ** 2 / (2 * (param[2] ** 2))) / (param[2] * np.sqrt(2 * np.pi)) + param[3] * x ** 2 + param[4] * x + \
        param[5]

def doFitPeak(xdata, ydata, odr = False, xerror = [], yerror = [], bkg = True, bkgForm = 'lin'):

    """
    Function for fitting single gaussian peak, with no or linear or quadratic background

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n
    bkg : boolean, optional
        indicates whether the background is included in the fit function\n
    bkgForm : str, optional
        indicates the type of background used, 'lin' for linear background, 'quad' for quadratic background\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            [
                'bk_a' : slope of linear background,
                'bk_b' : intercept of linear background,
            ], (if bkgForm == 'lin')
            [
                'bk_a' : quadratic term of quadratic background,
                'bk_b' : linear term of quadratic background,
                'bk_c' : constant of quadratic background,
            ], (if bkgForm == 'quad')
            'peak_amplitude': area of the gaussian peak,
            'peak_center': center of the gaussian peak,
            'peak_sigma': standard deviation of the gaussian peak,
            'peak_amplitude_err': error of area of the gaussian peak,
            'peak_center_err': error of center of the gaussian peak,
            'peak_sigma_err': error of standard deviation of the gaussian peak,
        }\n
    with fit function as\n
    `y =  peak_amplitude * exp(-(x - peak_center) ** 2 / (2 * peak_sigma ** 2)) / sqrt(2 * pi * peak_sigma ^ 2) + bk_a * x ** 2 + bk_b * x + bk_c`\n

    Raises
    ----------
    Exception
        when bkg is `True` the parameter bkgForm is neither `'lin'` nor `'quad'`
    """

    gModel = lmfit.models.GaussianModel(prefix = 'peak_')
    #Initial guess of gaussian peak parameters
    param1 = gModel.guess(ydata, x = xdata)
    bModel = None
    if bkg:
        ydatabkg = ydata - gaussianFunction([param1.valuesdict()['peak_amplitude'], param1.valuesdict()['peak_center'], \
            param1.valuesdict()['peak_sigma']], xdata)
        if bkgForm == 'lin':
            bModel = lmfit.models.LinearModel(prefix = 'bk_')
        elif bkgForm == 'quad':
            bModel = lmfit.models.QuadraticModel(prefix = 'bk_')
        else:
            raise Exception('doFitPeak: unable to parse background form \"' + bkgForm + '\"')
        #Initial guess of background parameters
        param2 = bModel.guess(ydatabkg, x = xdata)
        param = param1 + param2
    else:
        param = param1
    if odr:
        #ODR(Orthogonal distance regression) fit
        if len(xerror) == 0:
            if len(yerror) == 0:
                data = RealData(xdata, ydata)
            else:
                data = RealData(xdata, ydata, sy = yerror)
        else:
            if len(yerror) == 0:
                data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
            else:
                data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
        if bkg:
            if bkgForm == 'lin':
                model = Model(peakFunctionLin)
                paramInit = [param.valuesdict()['peak_amplitude'], param.valuesdict()['peak_center'], param.valuesdict()['peak_sigma'], \
                    param.valuesdict()['bk_slope'], param.valuesdict()['bk_intercept']]
            else:
                model = Model(peakFunctionQuad)
                paramInit = [param.valuesdict()['peak_amplitude'], param.valuesdict()['peak_center'], param.valuesdict()['peak_sigma'], \
                    param.valuesdict()['bk_a'], param.valuesdict()['bk_b'], param.valuesdict()['bk_c']]
        else:
            model = Model(gaussianFunction)
            paramInit = [param.valuesdict()['peak_amplitude'], param.valuesdict()['peak_center'], param.valuesdict()['peak_sigma']]
        odrFit = ODR(data, model, paramInit)
        odrFit.set_job(fit_type = 0)
        result = odrFit.run()
        fitResult = {
                        'peak_amplitude':           result.beta[0],
                        'peak_center':                result.beta[1],
                        'peak_sigma':                 result.beta[2],
                        'peak_amplitude_err':           result.sd_beta[0],
                        'peak_center_err':                result.sd_beta[1],
                        'peak_sigma_err':                 result.sd_beta[2],
            }
        if bkg:
            if bkgForm == 'lin':
                fitResult.update({'bk_a':       result.beta[3], 
                                          'bk_b':        result.beta[4], 
                    })
            else:
                fitResult.update({'bk_a':       result.beta[3], 
                                          'bk_b':        result.beta[4], 
                                          'bk_c':        result.beta[5], 
                    })
    else:
        #Least-squares fit
        if bkg:
            model = bModel + gModel
        else:
            model = gModel
        if len(yerror) == 0:
            result = model.fit(ydata, param, x = xdata)
        else:
            result = model.fit(ydata, param, x = xdata, weights = 1. / yerror)
        fitResult = {
                        'peak_amplitude':           result.best_values['peak_amplitude'],
                        'peak_center':                result.best_values['peak_center'],
                        'peak_sigma':                 result.best_values['peak_sigma'],
                        'peak_amplitude_err':           result.params['peak_amplitude'].stderr,
                        'peak_center_err':                result.params['peak_center'].stderr,
                        'peak_sigma_err':                 result.params['peak_sigma'].stderr,
            }
        if bkg:
            if bkgForm == 'lin':
                fitResult.update({'bk_a':       result.best_values['bk_slope'], 
                                          'bk_b':        result.best_values['bk_intercept'], 
                    })
            else:
                fitResult.update({'bk_a':       result.best_values['bk_a'], 
                                          'bk_b':        result.best_values['bk_b'], 
                                          'bk_c':        result.best_values['bk_c'], 
                    })
    return fitResult

def doFitDouble(xdata, ydata, fitRange, yerror = [], bkg = True, bkgForm = 'lin'):

    """
    CURRENTLY UNUSED AND NOT RECOMMENDED
    ----------

    Fit function for fitting overlapping double-gaussian peaks (namely for Co-60 data)\n
    ODR method is not used in this complex fit process, but can be added if necessary

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to be fit\n
    fitRange: list,
        containing the approximate ranges of the two peaks, in the form of `[peak1_lower, peak1_upper, peak2_lower, peak2_upper]`\
             in ascending order\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n
    bkg : boolean, optional
        indicates whether the background is included in the fit function\n
    bkgForm : str, optional
        indicates the type of background used, 'lin' for linear background, 'quad' for quadratic background\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            [
                'bk_a' : slope of linear background,
                'bk_b' : intercept of linear background,
            ], (if bkgForm == 'lin')
            [
                'bk_a' : quadratic term of quadratic background,
                'bk_b' : linear term of quadratic background,
                'bk_c' : constant of quadratic background,
            ], (if bkgForm == 'quad')
            'peak1_amplitude': area of the gaussian peak 1,
            'peak1_center': center of the gaussian peak 1,
            'peak1_sigma': standard deviation of the gaussian peak 1,
            'peak2_amplitude': area of the gaussian peak 2,
            'peak2_center': center of the gaussian peak 2,
            'peak2_sigma': standard deviation of the gaussian peak 2,
            'peak1_amplitude_err': error of area of the gaussian peak 1,
            'peak1_center_err': error of center of the gaussian peak 1,
            'peak1_sigma_err': error of standard deviation of the gaussian peak 1,
            'peak2_amplitude_err': error of area of the gaussian peak 2,
            'peak2_center_err': error of center of the gaussian peak 2,
            'peak2_sigma_err': error of standard deviation of the gaussian peak 2,
        }\n
    with fit function as\n
    `y = peak1_amplitude * exp(-(x - peak1_center) ^ 2 / (2 * peak1_sigma ^ 2)) / sqrt(2 * pi * peak1_sigma ^ 2) + \
        peak2_amplitude * exp(-(x - peak2_center) ^ 2 / (2 * peak2_sigma ^ 2)) / sqrt(2 * pi * peak2_sigma ^ 2)`\n

    Raises
    ----------
    Exception
        when bkg is `True` the parameter bkgForm is neither `'lin'` nor `'quad'`
    """

    gModel1 = lmfit.models.GaussianModel(prefix = 'peak1_')
    param1 = gModel1.guess(ydata[fitRange[0]:fitRange[1] + 1], x = xdata[fitRange[0]:fitRange[1] + 1])
    gModel2 = lmfit.models.GaussianModel(prefix = 'peak2_')
    param2 = gModel2.guess(ydata[fitRange[2]:fitRange[3] + 1], x = xdata[fitRange[2]:fitRange[3] + 1])
    bModel = None
    if bkg:
        ydatabkg = ydata - gaussianFunction([param1.valuesdict()['peak1_amplitude'], param1.valuesdict()['peak1_center'], \
            param1.valuesdict()['peak1_sigma']], xdata) - gaussianFunction([param2.valuesdict()['peak2_amplitude'], param2.valuesdict()\
            ['peak2_center'], param2.valuesdict()['peak2_sigma']], xdata)
        if bkgForm == 'lin':
            bModel = lmfit.models.LinearModel(prefix = 'bk_')
        elif bkgForm == 'quad':
            bModel = lmfit.models.QuadraticModel(prefix = 'bk_')
        else:
            raise Exception('doFitPeak: unable to parse background form \"' + bkgForm + '\"')
        #Initial guess of background parameters
        param3 = bModel.guess(ydatabkg, x = xdata)
        param = param1 + param2 + param3
        model = gModel1 + gModel2 + bModel
    else:
        param = param1 + param2
        model = gModel1 + gModel2
    if len(yerror) == 0:
        result = model.fit(ydata, param, x = xdata)
    else:
        result = model.fit(ydata, param, x = xdata, weights = 1. / yerror)
    fitResult = {
                    'peak1_amplitude':           result.best_values['peak1_amplitude'],
                    'peak1_center':                result.best_values['peak1_center'],
                    'peak1_sigma':                 result.best_values['peak1_sigma'],
                    'peak2_amplitude':           result.best_values['peak2_amplitude'],
                    'peak2_center':                result.best_values['peak2_center'],
                    'peak2_sigma':                 result.best_values['peak2_sigma'],
                    'peak1_amplitude_err':           result.params['peak1_amplitude'].stderr,
                    'peak1_center_err':                result.params['peak1_center'].stderr,
                    'peak1_sigma_err':                 result.params['peak1_sigma'].stderr,
                    'peak2_amplitude_err':           result.params['peak2_amplitude'].stderr,
                    'peak2_center_err':                result.params['peak2_center'].stderr,
                    'peak2_sigma_err':                 result.params['peak2_sigma'].stderr,
        }
    if bkg:
        if bkgForm == 'lin':
            fitResult.update({'bk_a':       result.best_values['bk_slope'], 
                                       'bk_b':       result.best_values['bk_intercept'], 
                })
        else:
            fitResult.update({'bk_a':       result.best_values['bk_a'], 
                                       'bk_b':       result.best_values['bk_b'], 
                                       'bk_c':       result.best_values['bk_c'], 
                })
    return fitResult

#-------------------------------------------------------Temperature-bias response fit functions-------------------------------------------------------

def quadFunction(param, x):

    """
    Auxiliary function to calculate the value of the quadratic function, used for ODR fit\n
    in the form of\n
    `f(x) = a * x ^ 2 + b * x + c`\n

    Parameters
    ----------
    param : list,
        parameters of the quadratic function, in the form of [a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the quadratic function\n
    """

    return param[0] * x * x + param[1] * x + param[2]

def doFitQuad(xdata, ydata, odr = False, xerror = [], yerror = []):

    """
    Function for fitting quadratic data, namely for temperature and bias response experiments

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': quadratic term,
            'fit_b': linear term,
            'fit_c': constant,
            'fit_a_err': error of quadratic term,
            'fit_b_err': error of linear term,
            'fit_c_err': error of constant,
        }
    with fit function as\n
    `y = fit_a * x ^ 2 + fit_b * x + fit_c`
    """

    qModel = lmfit.models.QuadraticModel(prefix = 'fit_')
    #Initial guess of the parameters
    param = qModel.guess(ydata, x = xdata)
    if odr:
        #ODR(Orthogonal distance regression) fit
        if len(xerror) == 0:
            if len(yerror) == 0:
               data = RealData(xdata, ydata)
            else:
               data = RealData(xdata, ydata, sy = yerror)
        else:
            if len(yerror) == 0:
                data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
            else:
                data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
        model = Model(quadFunction)
        odrFit = ODR(data, model, [param.valuesdict()['fit_a'], param.valuesdict()['fit_b'], param.valuesdict()['fit_c']])
        odrFit.set_job(fit_type = 0)
        result = odrFit.run()
        fitResult = {
                        'fit_a':            result.beta[0],
                        'fit_b':            result.beta[1],
                        'fit_c':            result.beta[2],
                        'fit_a_err':            result.sd_beta[0],
                        'fit_b_err':            result.sd_beta[1],
                        'fit_c_err':            result.sd_beta[2],
            }
    else:
        #Least-squares fit
        if len(yerror) == 0:
            result = qModel.fit(ydata, param, x = xdata)
        else:
            result = qModel.fit(ydata, param, x = xdata, weights = 1. / np.array(yerror))
        fitResult = {
                        'fit_a':            result.best_values['fit_a'],
                        'fit_b':            result.best_values['fit_b'],
                        'fit_c':            result.best_values['fit_c'],
                        'fit_a_err':            result.params['fit_a'].stderr,
                        'fit_b_err':            result.params['fit_b'].stderr,
                        'fit_c_err':            result.params['fit_c'].stderr,
            }
    return fitResult

def quad2DFunction(input, p00, p10, p01, p11, p20, p02, p21, p12, p22, p30, p31, p40):

    """
    Auxiliary function to calculate the value of 2D quadratic function for correlative temperature-bias fit\n
    in the form of\n
    `f(x, y) = p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ^ 2 + p02 * y ^ 2 + p21 * x ^ 2 * y + p12 * x * y ^ 2 + p22 * x ^ 2 * y ^ 2 + p30 * x ^ 3 + \
        p31 * x ^ 3 * y + p40 * x ^ 4`

    Parameters
    ----------
    input : array-like,
        array containing both x and y data, being `x = input[:, 0]` and `y = input[:, 1]`\n
    p00, p10, p01, p11, p20, p02, p21, p12, p22, p30, p31, p40 : float,
        corresponding parameters\n

    Returns
    ----------
    float or array-like,
        values of quadratic function corresponding to input x and y\n
    """

    xdata = input[:, 0]
    ydata = input[:, 1]
    return p00 + p10 * xdata + p01 * ydata + p11 * xdata * ydata + p20 * xdata ** 2 + p02 * ydata ** 2 + p21 * xdata ** 2 * ydata + \
        p12 *xdata * ydata ** 2 + p22 * xdata ** 2 * ydata ** 2 + p30 * xdata ** 3 + p31 * xdata ** 3 * ydata + p40 * xdata ** 4

def residualQuad2D(param, xdata, ydata, zdata):

    """
    Auxiliary function to calculate the residual of 2D quadratic surface model for correlative temperature-bias fit\n
    in the form of\n
    `res = zdata - f(xdata, ydata)`\n
    where\n
    `f(x, y) = p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ^ 2 + p02 * y ^ 2 + p21 * x ^ 2 * y + p12 * x * y ^ 2 + p22 * x ^ 2 * y ^ 2 + p30 * x ^ 3 + \
        p31 * x ^ 3 * y + p40 * x ^ 4`
    Parameters
    ----------
    param : list,
        parameters of 2D quadratic function, in the form of [p00, p10, p01, p11, p20, p02, p21, p12, p22, p30, p31, p40]\n
    xdata : float or array-like,
        value(s) of input x\n
    ydata : float or array-like,
        value(s) of input y\n
    zdata : float or array-like,
        value(s) of input z\n

    Returns
    ----------
    float or array-like,
        calculated residual
    """

    return zdata - quad2DFunction(np.dstack((xdata, ydata))[0], param[0], param[1], param[2], param[3], param[4], param[5], param[6], \
        param[7], param[8], param[9], param[10], param[11])

def doFitQuad2D(xdata, ydata, zdata, odr = False, xerror = [], yerror = [], zerror = []):

    """
    Function for fitting 2D quadratic data, namely for temperature-bias response experiments

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    zdata : array-like,
        the z data to fit\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n
    zerror : array-like, optional
        standard deviation of z data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'p00':                  constant term, 
            'p10':                  coefficient of term `x`, 
            'p01':                  coefficient of term `y`, 
            'p11':                  coefficient of term `x * y`, 
            'p20':                  coefficient of term `x ^ 2`, 
            'p02':                  coefficient of term `y ^ 2`, 
            'p21':                  coefficient of term `x ^ 2 * y`, 
            'p12':                  coefficient of term `x * y ^ 2`, 
            'p22':                  coefficient of term `x ^ 2 * y ^ 2`, 
            'p30':                  coefficient of term `x ^ 3`, 
            'p31':                  coefficient of term `x ^ 3 * y`, 
            'p40':                  coefficient of term `x ^ 4`, 
            'p00_err':            error of constant term, 
            'p10_err':            error of coefficient of term `x`, 
            'p01_err':            error of coefficient of term `y`, 
            'p11_err':            error of coefficient of term `x * y`, 
            'p20_err':            error of coefficient of term `x ^ 2`, 
            'p02_err':            error of coefficient of term `y ^ 2`, 
            'p21_err':            error of coefficient of term `x ^ 2 * y`, 
            'p12_err':            error of coefficient of term `x * y ^ 2`, 
            'p22_err':            error of coefficient of term `x ^ 2 * y ^ 2`, 
            'p30_err':            error of coefficient of term `x ^ 3`, 
            'p31_err':            error of coefficient of term `x ^ 3 * y`, 
            'p40_err':            error of coefficient of term `x ^ 4`, 
            'chisquare':         chi-square value to determine goodness of the fit, 
        }
    with fit function as\n
    `y = p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ^ 2 + p02 * y ^ 2 + p21 * x ^ 2 * y + p12 * x * y ^ 2 + p22 * x ^ 2 * y ^ 2 + p30 * x ^ 3 + \
        p31 * x ^ 3 * y + p40 * x ^ 4`
    """

    def quad2DFunctionInternal(param, x):

        """
        Auxiliary function to calculate the value of 2D quadratic function for correlative temperature-bias fit, used for ODR fit\n
        in the form of\n
        `f(x, y) = p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ^ 2 + p02 * y ^ 2 + p21 * x ^ 2 * y + p12 * x * y ^ 2 + p22 * x ^ 2 * y ^ 2 + p30 * x ^ 3 + \
            p31 * x ^ 3 * y + p40 * x ^ 4`

        Parameters
        ----------
        param : list, 
            parameters of the function, in the form of `[p00, p10, p01, p11, p20, p02, p21, p12, p22, p30, p31, p40]`\n
        x : array-like,
            array containing both x and y data, being `xdata = x[0]` and `ydata = x[1]`\n

        Returns
        ----------
        float or array-like,
            values of quadratic function corresponding to input x and y\n
        """

        xdata = x[0]
        ydata = x[1]
        p00 = param[0]
        p10 = param[1]
        p01 = param[2]
        p11 = param[3]
        p20 = param[4]
        p02 = param[5]
        p21 = param[6]
        p12 = param[7]
        p22 = param[8]
        p30 = param[9]
        p31 = param[10]
        p40 = param[11]
        return p00 + p10 * xdata + p01 * ydata + p11 * xdata * ydata + p20 * xdata ** 2 + p02 * ydata ** 2 + p21 * xdata ** 2 * ydata + \
            p12 *xdata * ydata ** 2 + p22 * xdata ** 2 * ydata ** 2 + p30 * xdata ** 3 + p31 * xdata ** 3 * ydata + p40 * xdata ** 4

    if odr:
        #ODR(Orthogonal distance regression) fit
        inputData = np.row_stack((xdata, ydata))
        if len(xerror) == 0:
            if len(zerror) == 0:
                data = RealData(inputData, zdata)
            else:
                data = RealData(inputData, zdata, sy = zerror)
        else:
            inputErr = np.row_stack((xerror, yerror))
            if len(zerror) == 0:
                data = RealData(inputData, zdata, sx = inputErr, fix = np.ones(inputErr.shape))
            else:
                data = RealData(inputData, zdata, sx = inputErr, sy = zerror, fix = np.ones(inputErr.shape))
        model = Model(quad2DFunctionInternal)
        odrFit = ODR(data, model, list(np.ones(12)))
        odrFit.set_job(fit_type = 0)
        result = odrFit.run()
        chisq = result.sum_square / (len(zdata) - len(result.beta))
        fitResult = {
                        'p00':            result.beta[0], 
                        'p10':            result.beta[1], 
                        'p01':            result.beta[2], 
                        'p11':            result.beta[3], 
                        'p20':            result.beta[4], 
                        'p02':            result.beta[5], 
                        'p21':            result.beta[6], 
                        'p12':            result.beta[7], 
                        'p22':            result.beta[8], 
                        'p30':            result.beta[9], 
                        'p31':            result.beta[10], 
                        'p40':            result.beta[11], 
                        'p00_err':            result.sd_beta[0], 
                        'p10_err':            result.sd_beta[1], 
                        'p01_err':            result.sd_beta[2], 
                        'p11_err':            result.sd_beta[3], 
                        'p20_err':            result.sd_beta[4], 
                        'p02_err':            result.sd_beta[5], 
                        'p21_err':            result.sd_beta[6], 
                        'p12_err':            result.sd_beta[7], 
                        'p22_err':            result.sd_beta[8], 
                        'p30_err':            result.sd_beta[9], 
                        'p31_err':            result.sd_beta[10], 
                        'p40_err':            result.sd_beta[11], 
                        'chisquare':          chisq, 
            }
    else:
        inputData = np.array(np.dstack((xdata, ydata))[0])
        popt, pcov = curve_fit(quad2DFunction, inputData, zdata, sigma = zerror, absolute_sigma = True, method = 'lm')
        params = popt.tolist()
        perr = np.sqrt(np.diag(pcov))
        chisq = sum((residualQuad2D(params, xdata, ydata, zdata) / zerror) ** 2) / (len(zdata) - len(params))
        fitResult = {
                        'p00':                   params[0], 
                        'p10':                   params[1], 
                        'p01':                   params[2], 
                        'p11':                   params[3], 
                        'p20':                   params[4], 
                        'p02':                   params[5], 
                        'p21':                   params[6], 
                        'p12':                   params[7],  
                        'p22':                   params[8], 
                        'p30':                   params[9], 
                        'p31':                   params[10], 
                        'p40':                   params[11], 
                        'p00_err':             perr[0], 
                        'p10_err':             perr[1], 
                        'p01_err':             perr[2], 
                        'p11_err':             perr[3], 
                        'p20_err':             perr[4], 
                        'p02_err':             perr[5], 
                        'p21_err':             perr[6], 
                        'p12_err':             perr[7], 
                        'p22_err':             perr[8], 
                        'p30_err':             perr[9], 
                        'p31_err':             perr[10], 
                        'p40_err':             perr[11], 
                        'chisquare':          chisq, 
            }
    return fitResult

#-------------------------------------------------------Leak current-temperature fit functions-------------------------------------------------------

#The functions used in this part is still not correct. It would be another TBD to re-write leak current fit part with codes from zxt

def expFunction(param, x):

    """
    Auxiliary function to calculate the value of the exponential function with a constant\n
    in the form of\n
    `f(x) = a * exp(x / b) + c`\n

    Parameters
    ----------
    param : list,
        parameters of the exponential function, in the form of [a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the exponential function
    """

    return param[0] * np.exp(x / param[1]) + param[2]

def doFitExp(xdata, ydata, odr = False, xerror = [], yerror = []):

    """
    Function for fitting exponential data for leak current-temperature

    Parameters
    ----------
    xdata : array-like or float,
        the x data to fit\n
    ydata : array-like or float,
        the y data to be fit\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like or float, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like or float, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude term,
            'fit_b': decay term,
            'fit_c': constant term,
            'fit_a_error': error of amplitude term,
            'fit_b_error': error of decay term,
            'fit_c_error': error of constant term,
        }
    with fit function as\n
    `y = fit_a * exp(x / fit_b) + fit_c`
    """

    eModel = lmfit.models.ExponentialModel(prefix = 'fit_')
    #Initial guess of the parameters
    param = eModel.guess(ydata, x = xdata)
    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, ydata)
        else:
            data = RealData(xdata, ydata, sy = yerror)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
    model = Model(expFunction)
    odrFit = ODR(data, model, [param.valuesdict()['fit_amplitude'], - param.valuesdict()['fit_decay'], 1e-6])
    if odr:
        #ODR(Orthogonal distance regression) fit
        odrFit.set_job(fit_type = 0)
    else:
        #Least-squares fit
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_c':            result.beta[2],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
                    'fit_c_err':            result.sd_beta[2],
        }
    return fitResult

def linExpFunction(param, x):

    """
    Auxiliary function for linear-exponential fitting of leak current-temperature curve, with some parameters fixed to theoretical values\n
    in the form of\n
    `f(x) = a * (b + T) * exp(- c * (1.1785 - 9.025e-5 * T - 3.05e-7 * T ^ 2) / T) + d`\n
    where `T = x + 273.15` represents absolute temperature

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [a, b, c, d]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the linear-exponential function
    """

    T = x + 273.15
    return param[0] * (param[1] + T) * np.exp(- param[2] * (1.1785 - 9.025e-5 * (T * np.ones(len(x))) - 3.05e-7 * (T * np.ones(len(x))) ** 2) / \
        (T * np.ones(len(x)))) + param[3]

def doLinExpFit(xdata, ydata, init, odr = False, xerror = [], yerror = []):

    """
    Function for linear-exponential fitting of leak current-temperature curve

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude,
            'fit_b': linear term,
            'fit_c': exponential term,
            'fit_d': constant term,
            'fit_a_err': error of amplitude,
            'fit_b_err': error of linear term,
            'fit_c_err': error of exponential term,
            'fit_d_err': error of constant term,
        }
    with fit function as\n
    `y = fit_a * (fit_b + T) * exp(- fit_c / T) + fit_d`\n
    where `T = x + 273.15` represents absolute temperature
    """

    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, ydata)
        else:
            data = RealData(xdata, ydata, sy = yerror)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
    model = Model(linExpFunction)
    odrFit = ODR(data, model, list(init))
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_c':            result.beta[2],
                    'fit_d':            result.beta[3],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
                    'fit_c_err':            result.sd_beta[2],
                    'fit_d_err':            result.sd_beta[3],
        }
    return fitResult

def revExpFunction(param, x):

    """
    Auxiliary function for exponential fitting with reversed exponential term of leak current-temperature curve\n
    in the form of\n
    `f(x) = a * exp(- b / T) + c`\n
    where `T = x + 273.15` represents absolute temperature

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the reversed exponential function
    """

    T = x + 273.15
    return param[0] * np.exp(- param[1] / (T * np.ones(len(x)))) + param[2]

def doRevExpFit(xdata, ydata, init, odr = False, xerror = [], yerror = []):

    """
    Function forexponential fitting with reversed exponential term of leak current-temperature curve

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude,
            'fit_b': exponential term,
            'fit_c': constant term,
            'fit_a_err': error of amplitude,
            'fit_b_err': error of exponential term,
            'fit_c_err': error of constant term,
        }
    with fit function as\n
    `y = fit_a * exp(- fit_b / T) + fit_c`\n
    where `T = x + 273.15` represents absolute temperature
    """

    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, ydata)
        else:
            data = RealData(xdata, ydata, sy = yerror)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
    model = Model(revExpFunction)
    odrFit = ODR(data, model, list(init))
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_c':            result.beta[2],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
                    'fit_c_err':            result.sd_beta[2],
        }
    return fitResult

def mixedExpFunction(param, x):

    """
    Auxiliary function for exponential fitting with mixed exponential term of leak current-temperature curve\n
    in the form of\n
    `f(x) = a * exp(- b * (1.1785 - 9.025e-5 * T - 3.05e-7 * T ^ 2) / T) + c`\n
    where `T = x + 273.15` represents absolute temperature

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the mixed exponential function
    """

    T = x + 273.15
    return param[0] * np.exp(- param[1] * (1.1785 - 9.025e-5 * T - 3.05e-7 * T ** 2) / T) + param[2]

def doMixedExpFit(xdata, ydata, init, odr = False, xerror = [], yerror = []):

    """
    Function forexponential fitting with mixed exponential term of leak current-temperature curve

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude
            'fit_b': exponential term
            'fit_c': constant term
            'fit_a_err': error of amplitude
            'fit_b_err': error of exponential term
            'fit_c_err': error of constant term
        }
    with fit function as\n
    `y = fit_a * exp(- fit_b * (1.1785 - 9.025e-5 * T - 3.05e-7 * T ^ 2) * / T) + fit_c`\n
    where `T = x + 273.15` represents absolute temperature
    """

    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, ydata)
        else:
            data = RealData(xdata, ydata, sy = yerror)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
    model = Model(mixedExpFunction)
    odrFit = ODR(data, model, list(init))
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_c':            result.beta[2],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
                    'fit_c_err':            result.sd_beta[2],
        }
    return fitResult

def fixedExpFunction(param, x):

    """
    Auxiliary function for fixed decay parameter exponential fitting of leak current-temperature curve\n
    in the form of\n
    `f(x) = a * exp(- b / T) + c`\n
    where `T = x + 273.15` represents absolute temperature, and \n
    `b = Eg(T)/2k = 1.1785 - 9.025e-5 * T - 3.05e-7 * T ^ 2` is set to fixed value

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the fixed-decay exponential function
    """

    T = x + 273.15
    return param[0] * np.exp(- 1.6e-19 * 1.1269 / (2 * 1.38e-23 * T)) + param[1]

def doFixedExpFit(xdata, ydata, init, odr = False, xerror = [], yerror = []):

    """
    Function for linear-exponential fitting of leak current-temperature curve

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude,
            'fit_b': constant term,
            'fit_a_err': error of amplitude,
            'fit_b_err': error of constant term,
        }
    with fit function as\n
    `y = fit_a * exp(- b / T) + fit_b`\n\n
    where `T = x + 273.15` represents absolute temperature\n
    and `b = Eg(T)/2k = 1.1785 - 9.025e-5 * T - 3.05e-7 * T ^ 2` ~ 1.1269eV @ 0-25 degrees celcius
    """

    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, ydata)
        else:
            data = RealData(xdata, ydata, sy = yerror)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
    model = Model(fixedExpFunction)
    odrFit = ODR(data, model, list(init))
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
        }
    return fitResult

#------------------------------------------------------------Leak current-bias fit functions------------------------------------------------------------

def revLinExpFunction(param, x):

    """
    Auxiliary function for exponential fitting with reversed exponential term and a linear term of leak current-bias curve\n
    in the form of\n
    `f(x) = a * (x - c) * exp(- b / (x - c)) + d`

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the reversed linear-exponential function
    """

    return param[0] * (x - param[2]) * np.exp(- param[1] / (x - param[2])) + param[3]

def doRevLinExpFit(xdata, ydata, init, odr = False, xerror = [], yerror = []):

    """
    Function forexponential fitting with reversed exponential term and a linear term of leak current-bias curve

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude,
            'fit_b': exponential term,
            'fit_c': linear term,
            'fit_d': constant term,
            'fit_a_err': error of amplitude,
            'fit_b_err': error of exponential term,
            'fit_c_err': error of linear term,
            'fit_d_err': error of constant term,
        }
    with fit function as\n
    `y = fit_a * (x - fit_c) * exp(- fit_b / (x - fit_c))`
    """

    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, ydata)
        else:
            data = RealData(xdata, ydata, sy = yerror)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
    model = Model(revLinExpFunction)
    odrFit = ODR(data, model, list(init))
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_c':            result.beta[2],
                    'fit_d':            result.beta[3],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
                    'fit_c_err':            result.sd_beta[2],
                    'fit_d_err':            result.sd_beta[3],
        }
    return fitResult

#--------------------------------------------------------------Energy-channel fit functions--------------------------------------------------------------

def linearFunction(param, x):

    """
    Auxiliary function to calculate the value of the linear function, used for ODR fit\n
    in the form of\n
    `f(x) = a * x + b`\n

    Parameters
    ----------
    param : list,
        parameters of the linear function, in the form of [a, b]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the linear function\n
    """

    return param[0] * x + param[1]

def doFitLin(xdata, ydata, odr = False, xerror = [], yerror = []):

    """
    Function for fitting linear data, namely for energy-channel fits

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': slope,
            'fit_b': intercept,
            'fit_a_err': error of slope,
            'fit_b_err': error of intercept,
        }
    with fit function as\n
    `y = fit_a * x + fit_b`
    """

    lModel = lmfit.models.LinearModel(prefix = 'fit_')
    #Initial guess of the parameters
    param = lModel.guess(ydata, x = xdata)
    if odr:
        #ODR(Orthogonal distance regression) fit
        if len(xerror) == 0:
            if len(yerror) == 0:
               data = RealData(xdata, ydata)
            else:
               data = RealData(xdata, ydata, sy = yerror)
        else:
            if len(yerror) == 0:
                data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
            else:
                data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
        model = Model(linearFunction)
        odrFit = ODR(data, model, [param.valuesdict()['fit_slope'], param.valuesdict()['fit_intercept']])
        odrFit.set_job(fit_type = 0)
        result = odrFit.run()
        fitResult = {
                        'fit_a':            result.beta[0],
                        'fit_b':            result.beta[1],
                        'fit_a_err':            result.sd_beta[0],
                        'fit_b_err':            result.sd_beta[1],
            }
    else:
        #Least-squares fit
        if len(yerror) == 0:
            result = lModel.fit(ydata, param, x = xdata)
        else:
            result = lModel.fit(ydata, param, x = xdata, weights = 1. / np.array(yerror))
        fitResult = {
                        'fit_a':            result.best_values['fit_slope'],
                        'fit_b':            result.best_values['fit_intercept'],
                        'fit_a_err':            result.params['fit_slope'].stderr,
                        'fit_b_err':            result.params['fit_intercept'].stderr,
            }
    return fitResult

def resolutionFunction(param, x):

    """
    THIS FUNCTION IS CURRENTLY UNUSED
    ----------

    Auxiliary function for energy resolution fit in EC part\n
    in the form of\n
    `f(x) = sqrt(a * x ^ 2 + b * x + c) / x` if a * x ^ 2 + b * x + c > 0\n
    `0` else

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [a, b, c]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    y : array-like or float,
        value(s) of the resolution fit function
    """

    inner = np.array(param[0] * x * x + param[1] * x + param[2])
    y = np.zeros(len(inner))
    qVal = np.where(inner > 0.)
    y[qVal] = np.sqrt(np.array(inner[qVal])) / np.array(x)[qVal]
    return y

def doResolutionFit(xdata, ydata, init, odr = False, xerror = [], yerror = []):

    """
    THIS FUNCTION IS CURRENTLY UNUSED
    ----------
    
    Function for energy resolution fit in EC part

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': quadratic term,
            'fit_b': linear term,
            'fit_c': constant,
            'fit_a_err': error of quadratic term,
            'fit_b_err': error of linear term,
            'fit_c_err': error of constant,
        }
    with fit function as\n
    `y = sqrt(fit_a * x ^ 2 + fit_b * x + fit_c) / x` if fit_a * x ^ 2 + fit_b * x + fit_c > 0\n
    `0` else
    """

    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, ydata)
        else:
            data = RealData(xdata, ydata, sy = yerror)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, ydata, sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, ydata, sx = xerror, sy = yerror, fix = np.ones(len(xerror)))
    model = Model(resolutionFunction)
    odrFit = ODR(data, model, list(init))
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_c':            result.beta[2],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
                    'fit_c_err':            result.sd_beta[2],
        }
    return fitResult

#------------------------------------------------------------Live time fit function------------------------------------------------------------

def convExpFunction(param, x):

    """
    Auxiliary function for live time fitting, the log value of distribution being reverse distribution of convolution of \
        43 exponential distributions\n
    in the form of\n
    `l(x) = ln(A * (x - 43 * C) ^ 42 * exp(- (x - 43 * C) / b) / b ^ 43)`\n
    ` = A + 42 * ln(x - 43 * C) - 43 * ln(b) - (x - 43 * C) / b`\n
    where `C = 50us` is the internal dead time of the MCU

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [A, b]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the convoluted exponential function
    """

    C = 50e-6
    return param[0] + 42.0 * np.log(x - 43 * C) - 43.0 * np.log(param[1]) - (x - 43 * C) / param[1]

def doConvExpFit(xdata, ydata, init, odr = False, xerror = [], yerror = []):

    """
    Function for convolution-of-exponential fitting of count rate correction

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude,
            'fit_b': deacy,
            'fit_a_err': error of amplitude,
            'fit_b_err': error of decay,
        }
    with fit function as\n
    `ln(y) = fit_a + 42 * ln(x - 43 * C) - 43 * ln(fit_b) - (x - 43 * C) / fit_b`\n
    where C = 50us is the internal dead time of the MCU
    """

    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, np.log(ydata))
        else:
            data = RealData(xdata, np.log(ydata), sy = yerror / ydata)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, np.log(ydata), sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, np.log(ydata), sx = xerror, sy = yerror / ydata, fix = np.ones(len(xerror)))
    model = Model(convExpFunction)
    odrFit = ODR(data, model, list(init))
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
        }
    return fitResult

def logLinFunction(param, x):

    """
    Auxiliary function for live time fitting with the form of linear function\n
    in the form of\n
    `l(x) = A - b * x`

    Parameters
    ----------
    param : list,
        parameters of the function, in the form of [A, b]\n
    x : array-like or float,
        input x value(s)\n

    Returns
    ----------
    array-like or float,
        value(s) of the linear function
    """

    return param[0] - param[1] * x

def doLogLinFit(xdata, ydata, odr = False, xerror = [], yerror = []):

    """
    Function for log-linear (exponential) fitting of count rate correction

    Parameters
    ----------
    xdata : array-like,
        the x data to fit\n
    ydata : array-like,
        the y data to fit\n
    init : list or array-like or float,
        the initial parameter values\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xerror : array-like, optional
        standard deviation of x data when using odr fit\n
    yerror : array-like, optional
        standard deviation of y data when using odr fit\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'fit_a': amplitude,
            'fit_b': average rate,
            'fit_a_err': error of amplitude,
            'fit_b_err': error of average rate,
        }
    with fit function as\n
    `ln(y) = fit_a - x * fit_b`
    """

    eModel = lmfit.models.ExponentialModel(prefix = 'fit_')
    param = eModel.guess(ydata, x = xdata)
    if len(xerror) == 0:
        if len(yerror) == 0:
            data = RealData(xdata, np.log(ydata))
        else:
            data = RealData(xdata, np.log(ydata), sy = yerror / ydata)
    else:
        if len(yerror) == 0:
            data = RealData(xdata, np.log(ydata), sx = xerror, fix = np.ones(len(xerror)))
        else:
            data = RealData(xdata, np.log(ydata), sx = xerror, sy = yerror / ydata, fix = np.ones(len(xerror)))
    model = Model(logLinFunction)
    odrFit = ODR(data, model, [np.log(param.valuesdict()['fit_amplitude']), 1. / param.valuesdict()['fit_decay']])
    if odr:
        odrFit.set_job(fit_type = 0)
    else:
        odrFit.set_job(fit_type = 2)
    result = odrFit.run()
    fitResult = {
                    'fit_a':            result.beta[0],
                    'fit_b':            result.beta[1],
                    'fit_a_err':            result.sd_beta[0],
                    'fit_b_err':            result.sd_beta[1],
        }
    return fitResult

#------------------------------------------------------------Spectrum and live time fit------------------------------------------------------------

def fitRateCorrect(filename, timeCorrect, plot = True, odr = False, rateStyle = ''):

    """
    Function for calculating correct total count rate and its error with the timeCorrect calculated when reading the data

    filename : str,
        name of the input data file\n
    timeCorrect : array-like or array of arrays,
        correct live time calculated when reading data\n
    plot : boolean, optional
        indicates the whether the figure of the fit result will be plotted\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    rateStyle : str, optional
        the style of calculating real count rate, '' for none, 's' for calculation with small data pack(512byte) live time data, 'p' for calculation \
            with large data packs(4096byte)\n

    Returns
    ----------
    rateAll : float,
        correct total count rate\n
    rateAllErr : float,
        error of correct total count rate\n

    Raises
    ----------
    Exception
        when rateStyle given is not in the available styles
    """

    styleAvailable = ['s', 'p'] #Option 'p' is currently NOT RECOMMENDED
    if not rateStyle in styleAvailable:
        raise Exception('fitRateCorrect: unknown count rate correction style ' + rateStyle)
    rateAll = 0.0
    rateAllErr = 0.0
    #Fill histogram of live time distribution
    specTime, xdata = np.histogram(timeCorrect, bins = 250, range=[np.min(timeCorrect), np.min([np.max(timeCorrect), parameters.liveTimeCut])])
    xdata = (xdata[:-1] + xdata[1:]) / 2
    C = parameters.deadTimeMCU
    center = 0.
    sigma = 0.
    if rateStyle == 's':
        #512-pack dead time correction
        qFit = (specTime > 0) * (xdata > C)
        result = doLogLinFit(xdata[qFit], specTime[qFit], odr = odr, yerror = gehrelsErr(specTime[qFit]))
        a = result['fit_a']
        b = result['fit_b']
        bErr = result['fit_b_err']
        rateAll = b
        rateAllErr = bErr
    else:
        #4096-pack dead time correction
        qFit = (specTime > 0) * (xdata > 43 * C)
        #Gaussian pre-fit for fit range correction
        res = doFitPeak(xdata[qFit], specTime[qFit], odr = False, bkg = False)
        center = res['peak_center']
        sigma = res['peak_sigma']
        #Fit range correction
        specTime, xdata = np.histogram(timeCorrect, bins = 100, range=[np.max([43 * C, center - 6. * sigma]), center + 6. * sigma])
        xdata = (xdata[:-1] + xdata[1:]) / 2
        qFit = (specTime > 0) * (xdata > 43 * C) * (xdata >= center - 4. * sigma) * (xdata <= center + 4. * sigma)
        result = doConvExpFit(xdata[qFit], specTime[qFit], (1., C), odr = odr, yerror = gehrelsErr(specTime[qFit]))
        a = result['fit_a']
        b = result['fit_b']
        bErr = result['fit_b_err']
        rateAll = 1 / b
        rateAllErr = bErr / b ** 2
    if plot:
        if rateStyle == 's':
            fitSpec = np.exp(logLinFunction([a, b], xdata))
            qPlot = xdata > C
        else:
            fitSpec = np.exp(convExpFunction([a, b], xdata))
            qPlot = xdata > 43 * C
        fig = plt.figure(figsize = (12, 8))
        gs = gridspec.GridSpec(1, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        ax = fig.add_subplot(gs[0])
        if rateStyle == 's':
            ax.step(xdata[qPlot], specTime[qPlot], where = 'mid', label = 'raw data', zorder = 1)
            ax.plot(xdata[qPlot], fitSpec[qPlot], label = 'Exponential Fit')
            ax.set_yscale('log')
            ax.set_ylim([1e-2, 1.2 * np.max(specTime)])
            ax.text(np.average(timeCorrect), max(specTime) * 1.1, 'live time = ' + str('%.2e' % (1 / b)) + ' $\pm$ ' + str('%.3e' % (bErr / b ** 2)) + \
                's\ncount rate = ' + str('%.2e' % (rateAll)) + ' $\pm$ ' + str('%.3e' % (rateAllErr)) + 'cps', fontsize = 10, bbox = dict(facecolor = 'pink', \
                alpha = 0.1), horizontalalignment = 'center', verticalalignment = 'center')
            ax.set_title('Fit of live time of ' + filename + '\nFit function: ' + r'$A e^{- b x}$')
        else:
            ax.step(xdata, specTime, where='mid', label = 'raw data', zorder = 1)
            ax.plot(xdata[qPlot], fitSpec[qPlot], label = 'Convoluted exponential Fit')
            ax.set_ylim([0, 1.2 * np.max(specTime)])
            ax.text(np.average(timeCorrect), max(specTime) * 1.1, 'live time = ' + str('%.2e' % (b)) + ' $\pm$ ' + str('%.3e' % (bErr)) + 's\ncount rate = ' + \
                str('%.2e' % (rateAll)) + ' $\pm$ ' + str('%.3e' % (rateAllErr)) + 'cps', fontsize = 10, bbox = dict(facecolor = 'pink', alpha = 0.1), \
                horizontalalignment = 'center', verticalalignment = 'center')
            ax.set_title('Fit of live time of ' + filename + '\nFit function: ' + r'$\frac{A (x - 43 C)^{42}}{(42)!b^{43}} e^{- (x - 43C) / b}$')
            #Vertiacl lines for marking
            ax.vlines(center - 3. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
            ax.vlines(center + 3. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
            ax.vlines(center - 4. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
            ax.vlines(center + 4. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
            ax.vlines(center - 5. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
            ax.vlines(center + 5. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
            ax.vlines(center - 6. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
            ax.vlines(center + 6. * sigma, 0, 1.2 * np.max(specTime), linewidth = 1, linestyle = ':')
        ax.set_xlabel('live time/s')
        ax.set_ylabel('count in bins')
        ax.legend(loc = 0)
        ax.grid()
        
        #Save/plot figure
        if parameters.saveFigNoPlot:
            filenameNoPath = extractFilename(filename)
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

#An idea to utilize error of correction factor: Use error propagation:
#   Let C be temp-bias correction factor, f(x) be original distribution, then
#   the corrected distribution is f1(x, C) = f(x / C)
#   Therofore the error propagation can be calculated as
#   sigma_f1 ** 2 = sigma_f ** 2 + (f'(x / C) * x / C ** 2) ** 2 * sigma_C ** 2
#Forget it. Won't work with the amplitude-calibration method currently used.
#Just do Monte-Carlo simulation to describe the error brought by temp-bias corretcion
def fitSpectrum(filename, amp, nbins, source, corr, time, fileOutput = False, singlech = False, bkg = False, odr = False, xRange = [],\
    channel = -1, corrErr = [], bkgAmp = [], bkgTime = -1.0, maxiter = 1, bound = 3.0, plot = True, rateStyle = '', rateAll = None, rateAllErr = None,\
    bkgRate = None, bkgRateErr = None, bkgForm = '', poissonErr = False, doCorr = True, lowerLim = -1.0, vLine = [], ecCorr = False, binWidth = None):

    """
    Function for fitting the spectrum

    Parameters
    ----------
    filename : str,
        name of input data file\n
    amp : array-like or array of arrays,
        ADC amplitude data\n
    nbins : int,
        number of bins to be used in the spectrum, within [1, parameters.adcMax]\n
    source : str,
        the source of the spectrum, currently supporting sources: `Am241`, `Ba133`, `Cs137`, `Na22`, `Co60`, `x`\n
    corr : list,
        temperature-bias correction factors for the data\n
    time : float,
        time taken to take the spectrum, usually calculated with event uscount\n
    fileOutput : boolean, optional
        the output style, False for figure output only, True for both figure and json file output\n
    singlech : boolean, optional
        indicates whether the data is single-channeled\n
    bkg : boolean, optional
        indicates whether the corresponding background is given. The corresponding total measurement times(source and background) should be \
            given(bkgTime) along with the background amplitude\n
    odr : boolean, optional
        indicates the fit method, if True then the fit is done with ODR method, or else the fit will be done with least-squares method\n
    xRange : list, optional
        specific fit range for x-ray data or uncalibrated data, in the form of `[lower, upper]`\n
    channel : int, optional
        channel number for single channel fits, within range [0, 3]\n
    corrErr : list, optional
        error of temperature-bias correction factors used for ODR fits, MAY be re-implemented for non-ODR fits\n
    bkgAmp : array-like or array of arrays, optional
        TEMP-BIAS CORRECTED amplitude of background data\n
    bkgTime : float, optional
        total measurement time of background data\n
    maxiter : int, optional
        maximum number of iterations, 0 for auto range correction off\n
    bound : float, optional
        boundary for auto rangecorrection in `\sigma`, with upper and lower bounds being `\mu - bound * \sigma` and `\mu + bound * \sigma`\n
    plot : boolean, optional
        indicates the whether the figure of the fit result will be plotted\n
    rateStyle : str, optional
        the style of calculating real count rate, '' for no correction, 's' for calculation with small data pack(512byte) live time data, 'p' for calculation \
            with large data packs(4096byte)\n
    rateAll : float or list, optional
        correct count rate of all spectrum in all 4 channels calculated when reading data, or count rates for 4 seperate channels, only used when \
            rateStyle is 's' or 'p'\n
    rateAllErr : float or list, optional
        error of rateAll\n
    bkgRate : float or list, optional
        correct count rate of all background spectrum in all 4 channels calculated when reading data, or count rates for 4 seperate channels, only \
            used when rateStyle is 's' or 'p'\n
    bkgRateErr : float or list, optional
        error of bkgRate\n
    bkgForm : str, optional
        indicates the type of internal peak background used, '' for no internal background, 'lin' for linear background, 'quad' for quadratic background\n
    poissonErr : boolean, optional
        indicates the method of calculating error of peak amplitude, True for poisson error, False for using error of fit model, when doing gaussian-only \
            fit\n
    doCorr : boolean, optional
        indicates whether the temperature-bias correction will be done, to avoid warning info in the output\n
    lowerLim : int, optional
        lower fit limit for fit range correction, specially for x-ray data fit\n
    vLine : list, optional
        positions for vertical lines to plot, specially for ec spectrum fits, and for SINGLE-CHANNELED DATA ONLY\n
    ecCorr : boolean, optional
        indicates whether the EC correction will be done\n
    binWidth : int or None, optional
        bin width to be used in the spectrum, within [1, parameters.adcMax], or None if not used\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
                'a':        amplitude,
                'b':        center,
                'c':        sigma,
                'a_err':    error of amplitude,
                'b_err':    error of center,
                'c_err':    error of sigma,
                'rate':        peak count rate,
                'rate_err':        error of peak count rate,
            }

    Raises
    ----------
    Exception
        when source given is not an available source, fit range is not given in correct format, channel out of bound when fitting single-channeled data, \
            background data not given when doing external background subtraction, internal background form is not an available form, error saving \
            output file or figure(s), error importing fit result from input file and error reading output config file

    NOTES
    ----------
    Suggested nbins for various sources at normal temperature and bias level: \n
    `Am241`, `Ba133`: 2048\n
    `Cs137`, `Na22`, `Co60`: 512\n
    `Th228`: 1024\n
    `x`: varies with different circumstances, try plotting raw data for a preliminary estimation\n
    MAY IMPLEMENT AN AUTOMATIC RANGE-FINDING WITH RESOLUTION FIT DATA\n
    EVALUATION OF ERROR INTRODUCED BY TEMP-BIAS CORRECTION MAY BE IMPLEMENTED, as a TBD
    """

    #Reference fit range and nbins data
    fitRange = parameters.fitRangeRef
    if filename in parameters.nbinsFileRef:
        nbins = parameters.nbinsFileRef[filename]
    if filename in parameters.binWidthRef:
        binWidth = parameters.binWidthRef[filename]

    #Read output config file
    config = {}
    try:
        with open(parameters.outputConfig, 'r') as fin:
            config = json.load(fin)
    except:
        raise Exception('fitSpectrum: unable to read file output configuration from file \"' + parameters.outputConfig + '\"')

    #For count rate correction, correct count rate calculated with fitRateCorrect should be specified
    styleAvailable = ['', 's', 'p']
    if not rateStyle in styleAvailable:
        raise Exception('fitSpectrum: unknown count rate correction style')
    elif rateAll is None or rateAllErr is None:
        raise Exception('fitSpectrum: corrected count rate or corresponding error not specified')
    if bkg and (bkgRate is None or bkgRateErr is None):
        raise Exception('fitSpectrum: background corrected count rate or corresponding error not specified')
    
    #Auto initial fit range correction using temp-bias correction factor
    rangeCorr = [1.0, 1.0, 1.0, 1.0]
    if not doCorr:
        rangeCorr = copy(corr)
        corr = [1.0, 1.0, 1.0, 1.0]

    #Source readout
    rangeLim = []
    if not source in fitRange:
        if source == 'x':
            #X-ray data, or specified initial fit range for specified source
            if singlech:
                if not (len(xRange) == 2 and xRange[0] >= 0 and xRange[1] > xRange[0]):
                    raise Exception('fitSpectrum: illegal format for fit range, please make sure the fit range is in the form of [lower, upper] and 0 <= lower '
                        '< upper')
                rangeLim = xRange
            else:
                if len(xRange) > 0:
                    for ich in range(4):
                        if not (len(xRange[ich]) == 2 and xRange[ich][0] >= 0 and xRange[ich][1] > xRange[ich][0]):
                            raise Exception('fitSpectrum: illegal format for fit range, please make sure the fit range is in the form of [lower, upper] and 0 <= lower '
                                '< upper')
                        rangeLim = xRange
        else:
            #Unknown source
            raise Exception('fitSpectrum: unknown source, supported sources: Am241, Ba133, Cs137, Na22, Co60, x')
    else:
        #Specified source
        if singlech:
            rangeLim = fitRange[source][channel]
            rangeLim[0] /= rangeCorr[channel]
            rangeLim[1] /= rangeCorr[channel]
            if len(xRange) > 0 and len(xRange[channel]) == 2:
                rangeLim = xRange[channel]
        else:
            rangeLim = fitRange[source]
            for ich in range(4):
                rangeLim[ich][0] /= rangeCorr[ich]
                rangeLim[ich][1] /= rangeCorr[ich]
            if len(xRange) > 0:
                for ich in range(4):
                    if len(xRange[ich]) == 2:
                        rangeLim[ich] = xRange[ich]

    #EC correction for fit range
    if ecCorr:
        rangeLim = ecCorrection(rangeLim, singlech, channel, dataCorr = False)

    #File output settings
    filenameNoPath = ''
    filenameNoAppend = ''
    saveFitPath = ''
    saveFigPath = ''
    if fileOutput and config['fitResult'] or plot and parameters.saveFigNoPlot or parameters.importFitResult:
        filenameNoPath = extractFilename(filename)
        filenameNoAppend = filenameNoPath[0:filenameNoPath.find(filenameNoPath.split('.')[-1]) - 1]
        if fileOutput or parameters.importFitResult:
            saveFitPath = parameters.saveFitPath + '/'
            if parameters.saveFitOwnPath:
                saveFitPath += (filenameNoAppend + '/')
            if parameters.importFitResult:
                fitResult = []
                try:
                    with open(saveFitPath + 'fit_' + filenameNoAppend + '.json', 'r') as fin:
                        fitResult = json.load(fin)
                except:
                    raise Exception('fitSpectrum: Error reading input file \"fit_' + filenameNoAppend + '.json\"')
                return fitResult
            else:
                if not os.path.isdir(saveFitPath):
                    os.makedirs(saveFitPath)
        if plot:
            saveFigPath = parameters.saveFigPath + '/'
            if parameters.saveFigOwnPath:
                saveFigPath += (filenameNoAppend + '/')
            if not os.path.isdir(saveFigPath):
                os.makedirs(saveFigPath)

    #Channel readout and spectrum creation with statistical error
    spectrumRaw = []
    specRange = []
    specWidth = None
    if singlech:
        if not isChannel(channel):
            raise Exception('fitSpectrum: channel number out of bound[0-3]')
        dataCorr = np.array(amp[channel]) * corr[channel]
        specRange = [0., parameters.adcMax * corr[channel]]
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech, channel)
            specRange = [0., parameters.energyPlotUpper]
        specWidth = specRange[1] - specRange[0]
        spectrumRaw, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth)
        spectrumStatErr = gehrelsErr(spectrumRaw)
    else:
        specRange = []
        specWidth = []
        for ich in range(4):
            amp[ich] = np.array(amp[ich])
            specRange.append([0., parameters.adcMax * corr[ich]])
            specWidth.append(specRange[ich][1] - specRange[ich][0])
        dataCorr = amp * np.array(corr)
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech)
            specRange = [0., parameters.energyPlotUpper]
            specWidth = specRange[1] - specRange[0]
        spectrumRaw, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth)
        spectrumStatErr = []
        for ich in range(4):
            spectrumStatErr.append(gehrelsErr(spectrumRaw[ich]))

    #Rate correction
    rateFactor = 1.0
    rateFactorErr = 0.0
    if rateStyle != '':
        if isinstance(rateAll, list):
            rateFactor = []
            rateFactorErr = []
            for ich in range(4):
                rateFactor.append(rateAll[ich] / float(len(amp[ich])))
                rateFactorErr.append(rateAllErr[ich] / float(len(amp[ich])))
        else:
            countAll = 0.0
            for ich in range(4):
                countAll += float(len(amp[ich]))
            rateFactor = rateAll / countAll
            rateFactorErr = rateAllErr / countAll
        if singlech:
            if isinstance(rateAll, list):
                spectrumErr = np.sqrt((spectrumRaw * rateFactorErr[channel]) ** 2 + (spectrumStatErr * rateFactor[channel]) ** 2)
                spectrumRaw = np.array(spectrumRaw) * rateFactor[channel]
            else:
                spectrumErr = np.sqrt((spectrumRaw * rateFactorErr) ** 2 + (spectrumStatErr * rateFactor) ** 2)
                spectrumRaw = np.array(spectrumRaw) * rateFactor
        else:
            if isinstance(rateAll, list):
                spectrumErr = []
                for ich in range(4):
                    spectrumErr.append([])
                    spectrumErr[ich] = np.sqrt((spectrumRaw[ich] * rateFactorErr[ich]) ** 2 + (spectrumStatErr[ich] * rateFactor[ich]) ** 2)
                spectrumErr = np.array(spectrumErr)
                for ich in range(4):
                    spectrumRaw[ich] = np.array(spectrumRaw[ich]) * rateFactor[ich]
            else:
                spectrumErr = []
                for ich in range(4):
                    spectrumErr.append([])
                    spectrumErr[ich] = np.sqrt((spectrumRaw[ich] * rateFactorErr) ** 2 + (spectrumStatErr[ich] * rateFactor) ** 2)
                spectrumErr = np.array(spectrumErr)
                for ich in range(4):
                    spectrumRaw[ich] = np.array(spectrumRaw[ich]) * rateFactor
    else:
        spectrumRaw = np.array(spectrumRaw) / time
        spectrumErr = np.array(spectrumStatErr) / time

    #External background data spectrum and correction
    bkgRateFactor = 1.
    bkgRateFactorErr = 0.
    spectrum = []
    bkgSpectrum = []
    if bkg:
        if len(bkgAmp) == 0 or rateStyle == '' and bkgTime < 0:
            raise Exception('fitSpectrum: background data not given')
        if rateStyle == '':
            bkgRateFactor = 1. / bkgTime
        else:
            if isinstance(bkgRate, list):
                bkgRateFactor = []
                bkgRateFactorErr = []
                for ich in range(4):
                    bkgRateFactor.append(bkgRate[ich] / float(len(bkgAmp[ich])))
                    bkgRateFactorErr.append(bkgRateErr[ich] / float(len(bkgAmp[ich])))
            else:
                bkgCountAll = 0.0
                for ich in range(4):
                    bkgCountAll += float(len(bkgAmp[ich]))
                bkgRateFactor = bkgRate / bkgCountAll
                bkgRateFactorErr = bkgRateErr / bkgCountAll
        if singlech:
            bkgSpectrum, _ = getSpectrum(bkgAmp[channel], nbins, singlech, specRange, binWidth = binWidth)
            bkgSpectrumStatErr = gehrelsErr(bkgSpectrum)
            #Subtraction of external background
            if isinstance(bkgRate, list):
                spectrum = spectrumRaw - bkgSpectrum * bkgRateFactor[channel]
                spectrumErr = np.sqrt(spectrumErr ** 2 + (bkgSpectrum * bkgRateFactorErr[channel]) ** 2 + (bkgSpectrumStatErr * bkgRateFactor[channel]) ** 2)
                bkgSpectrum = np.array(bkgSpectrum) * bkgRateFactor[channel]
            else:
                spectrum = spectrumRaw - bkgSpectrum * bkgRateFactor
                spectrumErr = np.sqrt(spectrumErr ** 2 + (bkgSpectrum * bkgRateFactorErr) ** 2 + (bkgSpectrumStatErr * bkgRateFactor) ** 2)
                bkgSpectrum = np.array(bkgSpectrum) * bkgRateFactor
        else:
            specRange = [0., parameters.adcMax * corr[channel]]
            bkgSpectrum, _ = getSpectrum(bkgAmp, nbins, singlech, specRange, binWidth = binWidth)
            bkgSpectrumStatErr = []
            for ich in range(4):
                bkgSpectrumStatErr.append(gehrelsErr(bkgSpectrum[ich]))
            #Subtraction of external background
            if isinstance(bkgRate, list):
                for ich in range(4):
                    spectrum.append(spectrumRaw[ich] - bkgSpectrum[ich] * bkgRateFactor[ich])
                    spectrumErr[ich] = np.sqrt(spectrumErr[ich] ** 2 + (bkgSpectrum[ich] * bkgRateFactorErr[ich]) ** 2 + (bkgSpectrumStatErr[ich] * \
                        bkgRateFactor[ich]) ** 2)
                    bkgSpectrum[ich] = np.array(bkgSpectrum[ich]) * bkgRateFactor[ich]
            else:
                for ich in range(4):
                    spectrum.append(spectrumRaw[ich] - bkgSpectrum[ich] * bkgRateFactor)
                    spectrumErr[ich] = np.sqrt(spectrumErr[ich] ** 2 + (bkgSpectrum[ich] * bkgRateFactorErr) ** 2 + (bkgSpectrumStatErr[ich] * bkgRateFactor) \
                        ** 2)
                    bkgSpectrum[ich] = np.array(bkgSpectrum[ich]) * bkgRateFactor
            spectrum = np.array(spectrum)
    else:
        spectrum = spectrumRaw

    #Readout of interal background type
    bkgFormAvailable = ['', 'lin', 'quad']
    if not bkgForm in bkgFormAvailable:
        raise Exception('fitSpectrum: unknown internal background type \"' + bkgForm + '\"')

    a = 0.
    b = 0.
    c = 0.
    minDiff = 1e-4
    #single channel fits
    if singlech:
        ich = channel
        if ich < 0 or ich > 3:
            raise Exception('fitSpectrum: channel number out of bound[0-3]')
        #Divide bin counts by bin width
        if binWidth is None:
            spectrum = spectrum * nbins / specWidth
            spectrumErr = spectrumErr * nbins / specWidth
        else:
            spectrum = spectrum / binWidth
            spectrumErr = spectrumErr / binWidth
        q1 = (x >= rangeLim[0]) * (x <= rangeLim[1])
        result = doFitPeak(x[q1], spectrum[q1], odr, yerror = spectrumErr[q1], bkg = (bkgForm != ''), bkgForm = bkgForm)
        amplitude = result['peak_amplitude']
        center = result['peak_center']
        sigma = result['peak_sigma']
        if bkgForm == 'lin':
            a = result['bk_a']
            b = result['bk_b']
        elif bkgForm == 'quad':
            a = result['bk_a']
            b = result['bk_b']
            c = result['bk_c']
        lastCenter, lastSigma = center, sigma

        #Iteration part
        for iit in range(maxiter):
            q2 = (x >= max(center - bound * sigma, lowerLim)) * (x <= center + bound * sigma)
            result = doFitPeak(x[q2], spectrum[q2], odr, yerror = spectrumErr[q2], bkg = (bkgForm != ''), bkgForm = bkgForm)
            amplitude = result['peak_amplitude']
            center = result['peak_center']
            sigma = result['peak_sigma']
            if bkgForm == 'lin':
                a = result['bk_a']
                b = result['bk_b']
            elif bkgForm == 'quad':
                a = result['bk_a']
                b = result['bk_b']
                c = result['bk_c']
            if np.abs((lastCenter - center) / lastCenter) < minDiff and np.abs((lastSigma - sigma) / lastSigma) < minDiff:
                break
            lastCenter, lastSigma = center, sigma

        ampErr = result['peak_amplitude_err']
        centerErr = result['peak_center_err']
        sigmaErr = result['peak_sigma_err']
        resolution = 2 * np.sqrt(2 * np.log(2)) * sigma / center
        resolutionErr = 2 * np.sqrt(2 * np.log(2)) * np.sqrt((sigmaErr / center) ** 2 + (sigma * centerErr / center ** 2) ** 2)
        #Count rate calculation part
        if rateStyle == '':
            rate = amplitude
            if bkgForm != '' or poissonErr:
                rateErr = np.sqrt(amplitude / time)
            else:
                rateErr = ampErr
        else:
            rate = amplitude
            if bkgForm != '' or poissonErr:
                if isinstance(rateAll, list):
                    rateErr = np.sqrt(amplitude * rateFactor[ich] + (amplitude * rateFactorErr[ich] / rateFactor[ich]) ** 2)
                else:
                    rateErr = np.sqrt(amplitude * rateFactor + (amplitude * rateFactorErr / rateFactor) ** 2)
            else:
                rateErr = ampErr

        #Fit result
        fitResult = {
                                        'a':            amplitude,
                                        'b':            center,
                                        'c':            sigma,
                                        'a_err':        ampErr,
                                        'b_err':        centerErr,
                                        'c_err':        sigmaErr,
                                        'rate':        rate,
                                        'rate_err':        rateErr,
                }

        #Plot part
        if plot:
            fitPeak = gaussianFunction([amplitude, center, sigma], x)
            fitTotal = copy(fitPeak)
            fitBk = None
            if bkgForm == 'lin':
                fitBk = a * x + b
                fitTotal += fitBk
            elif bkgForm == 'quad':
                fitBk = quadFunction([a, b, c], x)
                fitTotal += fitBk
            #Plot range
            qplot = (x >= max(center - bound * sigma, lowerLim)) * (x < center + bound * sigma)
            fig = plt.figure(figsize = (12, 8))
            gs = gridspec.GridSpec(2, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95, height_ratios = [4, 1])
            ax = fig.add_subplot(gs[0])
            ax1 = fig.add_subplot(gs[1])
            if bkg:
                if binWidth is None:
                    spectrumRaw *= nbins / specWidth
                    bkgSpectrum *= nbins / specWidth
                else:
                    spectrumRaw /= binWidth
                    bkgSpectrum /= binWidth
                ax.step(x, spectrumRaw, where = 'mid', label = 'raw data', zorder = 1)
                ax.step(x, bkgSpectrum, where = 'mid', label = 'background', zorder = 1)
                ax.step(x, spectrum, where = 'mid', label = 'data without bkg', zorder = 1)
            else:
                ax.step(x, spectrum, where = 'mid', label = 'raw data', zorder = 1)
            ax.plot(x[qplot], fitPeak[qplot], label = 'Gaussian Peak Fit')
            if bkgForm != '':
                ax.plot(x[qplot], fitBk[qplot], label = 'Background Fit')
                ax.plot(x[qplot], fitTotal[qplot], label = 'Gaussian Peak and Background Fit')
            for iline in range(len(vLine)):
                ax.vlines(vLine[iline], 0, np.max(spectrum), linewidth = 1, linestyle = ':')
            residual = fitTotal[qplot] - spectrum[qplot]
            ax.text(center - 0.5 * sigma, max(fitPeak[qplot]) * 0.25, 'count rate = ' + str('%.2e' % (rate)) + ' $\pm$ ' + str('%.3e' % (rateErr)) + '\ncenter = ' + \
                str('%.3e' % center) + ' $\pm$ ' + str('%.3e' % centerErr) + '\nresolution = ' + str('%.3e' % resolution) + ' $\pm$ ' + str('%.3e' % resolutionErr), \
                fontsize = 10, bbox = dict(facecolor = 'pink', alpha = 0.1), horizontalalignment = 'center', verticalalignment = 'center')
            ax.set_ylim([0, 1.2 * np.max(spectrum[qplot])])
            ax1.plot(x[qplot], residual)
            #Lower and upper bounds for plot
            plotLowerLim = 0.
            if parameters.plotSpecLogScale:
                plotLowerLim = parameters.logScaleBegin
            lower = center - 6 * sigma
            if lower < plotLowerLim:
                lower = plotLowerLim
            upper = center + 6 * sigma
            if upper > parameters.adcMax * corr[ich]:
                upper = parameters.adcMax * corr[ich]
            ax.set_xlim([lower, upper])
            if parameters.plotSpecLogScale:
                ax.set_xscale('log')
                ax.set_yscale('log')
            ax.set_title('Spectrum fit of data from channel ' + str(ich) + ' in ' + filename + '\nFit function: ' + \
                r'$\frac{a}{\sqrt{2\pi}\sigma} e^{[{-{(x - \mu)^2}/{{2\sigma}^2}}]}$')
            ax.set_xlabel('ADC/channel')
            ax.set_ylabel('count rate/cps')
            ax1.set_ylabel('residual/cps')
            ax.legend(loc = 0)
            ax.grid()
            ax1.grid()
            if parameters.saveFigNoPlot:
                try:
                    fig.savefig('{}/{}_spectrum_fit_ch{}.png'.format(saveFigPath, filenameNoAppend, channel), transparent = True)
                    plt.close(fig)
                except:
                    print('fitRateCorrect: error saving spectrum fit figure for file \"' + filenameNoPath + '\"')
            else:
                plt.show()

        #File output
        if fileOutput and config['fitResult']:
            try:
                with open(saveFitPath + 'fit_' + filenameNoAppend + '_ch' + str(channel) + '.json', 'w') as fout:
                    json.dump(fitResult, fout)
            except:
                raise Exception('fitSpectrum: Error writing output file \"fit_' + filenameNoAppend + '.json\"')

    #Multiple channel fits
    else:
        fig = None
        gs = None
        if plot:
            fig = plt.figure(figsize = (12, 8))
            gs = gridspec.GridSpec(4, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        fitResult = []

        for ich in range(4):
            #Divide bin counts by bin width
            if binWidth is None:
                if isinstance(specWidth, list):
                    spectrum[ich] = spectrum[ich] * nbins / specWidth[ich]
                    spectrumErr[ich] = spectrumErr[ich] * nbins / specWidth[ich]
                else:
                    spectrum[ich] = spectrum[ich] * nbins / specWidth
                    spectrumErr[ich] = spectrumErr[ich] * nbins / specWidth
            else:
                spectrum[ich] = spectrum[ich] / binWidth
                spectrumErr[ich] = spectrumErr[ich] / binWidth
            q1 = (x[ich] >= rangeLim[ich][0]) * (x[ich] <= rangeLim[ich][1]) #Initial ROI
            result = doFitPeak(x[ich][q1], spectrum[ich][q1], odr, yerror = spectrumErr[ich][q1], bkg = (bkgForm != ''), bkgForm = bkgForm)
            amplitude = result['peak_amplitude']
            center = result['peak_center']
            sigma = result['peak_sigma']
            if bkgForm == 'lin':
                a = result['bk_a']
                b = result['bk_b']
            elif bkgForm == 'quad':
                a = result['bk_a']
                b = result['bk_b']
                c = result['bk_c']
            lastCenter, lastSigma = center, sigma

            #Iteration part
            for iit in range(maxiter):
                q2 = (x[ich] >= max(center - bound * sigma, lowerLim)) * (x[ich] <= center + bound * sigma)
                result = doFitPeak(x[ich][q2], spectrum[ich][q2], odr, yerror = spectrumErr[ich][q2], bkg = (bkgForm != ''), bkgForm = bkgForm)
                amplitude = result['peak_amplitude']
                center = result['peak_center']
                sigma = result['peak_sigma']
                if bkgForm == 'lin':
                    a = result['bk_a']
                    b = result['bk_b']
                elif bkgForm == 'quad':
                    a = result['bk_a']
                    b = result['bk_b']
                    c = result['bk_c']
                if np.abs((lastCenter - center) / lastCenter) < minDiff and np.abs((lastSigma - sigma) / lastSigma) < minDiff:
                    break
                lastCenter, lastSigma = center, sigma

            ampErr = result['peak_amplitude_err']
            centerErr = result['peak_center_err']
            sigmaErr = result['peak_sigma_err']
            resolution = 2 * np.sqrt(2 * np.log(2)) * sigma / center
            resolutionErr = 2 * np.sqrt(2 * np.log(2)) * np.sqrt((sigmaErr / center) ** 2 + (sigma * centerErr / center ** 2) ** 2)
            #Count rate calculation part
            if rateStyle == '':
                rate = amplitude
                if bkgForm != '' or poissonErr:
                    rateErr = np.sqrt(amplitude / time)
                else:
                    rateErr = ampErr
            else:
                rate = amplitude
                if bkgForm != '' or poissonErr:
                    if isinstance(rateAll, list):
                        rateErr = np.sqrt(amplitude * rateFactor[ich] + (amplitude * rateFactorErr[ich] / rateFactor[ich]) ** 2)
                    else:
                        rateErr = np.sqrt(amplitude * rateFactor + (amplitude * rateFactorErr / rateFactor) ** 2)
                else:
                    rateErr = ampErr

            #Fit result
            fitResult.append({
                                        'a':            amplitude,
                                        'b':            center,
                                        'c':            sigma,
                                        'a_err':        ampErr,
                                        'b_err':        centerErr,
                                        'c_err':        sigmaErr,
                                        'rate':        rate,
                                        'rate_err':        rateErr,
                })

            #Plot part
            if plot:
                fitPeak = gaussianFunction([amplitude, center, sigma], x[ich])
                fitTotal = copy(fitPeak)
                fitBk = None
                if bkgForm == 'lin':
                    fitBk = a * x[ich] + b
                    fitTotal += fitBk
                elif bkgForm == 'quad':
                    fitBk = quadFunction([a, b, c], x[ich])
                    fitTotal += fitBk
                #Plot range
                qplot = (x[ich] >= max(center - bound * sigma, lowerLim)) * (x[ich] < center + bound * sigma)
                ax = fig.add_subplot(gs[ich])
                if bkg:
                    if binWidth is None:
                        if isinstance(specWidth, list):
                            spectrumRaw[ich] *= nbins / specWidth[ich]
                            bkgSpectrum[ich] *= nbins / specWidth[ich]
                        else:
                            spectrumRaw[ich] *= nbins / specWidth
                            bkgSpectrum[ich] *= nbins / specWidth
                    else:
                        spectrumRaw[ich] /= binWidth
                        bkgSpectrum[ich] /= binWidth
                    plt.step(x[ich], spectrumRaw[ich], where = 'mid', label = 'raw data', zorder = 1)
                    plt.step(x[ich], bkgSpectrum[ich], where = 'mid', label = 'background', zorder = 1)
                    plt.step(x[ich], spectrum[ich], where = 'mid', label = 'data without bkg', zorder = 1)
                else:
                    plt.step(x[ich], spectrum[ich], where = 'mid', label = 'raw data', zorder = 1)
                plt.plot(x[ich][qplot], fitPeak[qplot], label = 'Gaussian Peak Fit')
                if bkgForm != '':
                    plt.plot(x[ich][qplot], fitBk[qplot], label = 'Background Fit')
                    plt.plot(x[ich][qplot], fitTotal[qplot], label = 'Gaussian Peak and Background Fit')
                plt.text(center - 0.5 * sigma, max(fitPeak[qplot]) * 0.25, 'count rate = ' + str('%.2e' % (rate)) + ' $\pm$ ' + str('%.3e' % (rateErr)) + '\ncenter = ' + \
                    str('%.3e' % center) + ' $\pm$ ' + str('%.3e' % centerErr) + '\nresolution = ' + str('%.3e' % resolution) + ' $\pm$ ' + str('%.3e' % resolutionErr), \
                    fontsize = 10, bbox = dict(facecolor = 'pink', alpha = 0.1), horizontalalignment = 'center', verticalalignment = 'center')
                ax.set_ylim([0, 1.2 * np.max(spectrum[ich][qplot])])
                #Lower and upper bounds for plot
                plotLowerLim = 0.
                if parameters.plotSpecLogScale:
                    plotLowerLim = parameters.logScaleBegin
                lower = center - 6 * sigma
                if lower < plotLowerLim:
                    lower = plotLowerLim
                upper = center + 6 * sigma
                if upper > parameters.adcMax * corr[ich]:
                    upper = parameters.adcMax * corr[ich]
                ax.set_xlim([lower, upper])
                if parameters.plotSpecLogScale:
                    ax.set_xscale('log')
                    ax.set_yscale('log')
                if ich == 0:
                    ax.set_title('Spectrum fit of data from ' + filename + '\nFit function: ' + r'$\frac{a}{\sqrt{2\pi}\sigma} e^{[{-{(x - \mu)^2}/{{2\sigma}^2}}]}$')
                ax.set_xlabel('ADC/channel')
                ax.set_ylabel('count rante/cps')
                ax.legend(loc = 0)
                ax.grid()
        if plot:
            if parameters.saveFigNoPlot:
                try:
                    fig.savefig('{}/{}_spectrum_fit.png'.format(saveFigPath, filenameNoAppend), transparent = True)
                    plt.close(fig)
                except:
                    print('fitRateCorrect: error saving spectrum fit figure for file \"' + filenameNoPath + '\"')
            else:
                plt.show()

        #File output
        if fileOutput and config['fitResult']:
            try:
                with open(saveFitPath + 'fit_' + filenameNoAppend + '.json', 'w') as fout:
                    json.dump(fitResult, fout)
            except:
                raise Exception('fitSpectrum: Error writing output file \"fit_' + filenameNoAppend + '.json\"')

    return fitResult
