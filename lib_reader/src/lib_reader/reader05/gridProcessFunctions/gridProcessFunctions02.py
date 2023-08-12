# -*- coding:utf-8 -*-

"""
GRID02 I/O and Basic Process Functions Library
----------
Basic data readout/output and process functions\n
`v1.0.0` GRID02 and GRID02-based detector calibration result analysis, also quick view of GRID02 in-orbit data, by ghz, cjr and ydx\n
"""

#**************************************************************************************************************************************************
#***************************************************************Import packages***************************************************************
#**************************************************************************************************************************************************

import gridParametersCommon as parameters
import gridParameters.gridParameters02 as specificParam
import gridBasicFunctions as basic
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import struct
import os
import re
import json, csv

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
        name of the config file (with full path), namely `xrayConfig`, `sourceConfig` and `HPGeConfig` in gridParameters02.py\n
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
        name of the config file (with full path), namely `xrayProcessConfig` and `sourceProcessConfig` in gridParameters02.py\n
    
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

#****************************************************************************************************************************************************************
#*******************************************************************Data readout functions*****************************************************************
#****************************************************************************************************************************************************************

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
    filenameNoPath = basic.extractFilename(filename)

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

    print('dataReadout: processing \"' + filenameNoPath + '\"')

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
        elif filenameNoPath in specificParam.reqFileRef:
            telreq = specificParam.reqFileRef[filenameNoPath][0]
            scireq = specificParam.reqFileRef[filenameNoPath][1]
        else:
            telreq = -1
            scireq = -1


    #------------------------------------------------------------Readout for binary files------------------------------------------------------------
    if isBin:
        #Patterns of science and telemetry data packs
        sciPattern = re.compile(specificParam.pattrenRef['binSci'], re.S)
        telPattern = re.compile(specificParam.pattrenRef['binTel'], re.S)
        #Find positions of all science and telemetry data packs
        sciPackPos = basic.findPackPos(lines, sciPattern)
        telPackPos = basic.findPackPos(lines, telPattern)
        #Get all science and telemetry data packs
        sciPackData = []
        telPackData = []
        for idata in range(len(sciPackPos)):
            sciPackData.append(lines[sciPackPos[idata]:sciPackPos[idata] + 512])
        for idata in range(len(telPackPos)):
            telPackData.append(lines[telPackPos[idata]:telPackPos[idata] + 512])
        #Do deduplicate data and generate deduplicated raw data, without readout
        if specificParam.deduplicateFile:
            lines = basic.deduplicate(sciPackData, telPackData, specificParam.deduplicateFile)
            filenameNoAppend = filenameNoPath[0:filenameNoPath.find(filenameNoPath.split('.')[-1]) - 1]
            with open(filenameNoAppend + '_deduplicated.' + filenameNoPath.split('.')[-1], 'wb+') as fout:
                fout.write(lines)
        else:
            if doDeduplicate:
                #Deduplicate data
                sciPackData, telPackData = basic.deduplicate(sciPackData, telPackData)

            utcLaunch = 0.
            if parameters.utcCheck:
                utcLaunch = basic.convertStrToTimestamp(specificParam.launchTime)
            telcuruscount = 0
            #Readout of telemetry data packs
            for telData in telPackData:
                if len(telData) < 498:
                    #Avoid error caused by incomplete data packs
                    crcError += 1#49140 wrt this line, 53173 with this line, actual 49341
                    continue
                #CRC check
                if not basic.crcCheck(telData[:496], telData[496:498]):
                    crcError += 1
                    continue
                for it in range(7):
                    #Get current uscount value
                    telcuruscount = struct.unpack('>Q', telData[15 + 70 * it:23 + 70 * it])[0] / specificParam.internalFreq
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
                            iMon[ich].append(float(struct.unpack('>H', telData[47 + 2 * ich + 70 * it:49 + 2 * ich + 70 * it])[0]) / 4096.0 * 3.3 / specificParam.internalResistance)
                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * specificParam.internalResistance)
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
                if not basic.crcCheck(sciData[:510], sciData[510:512]):
                    crcError += 1
                    continue
                #Channel of currrent event
                ch = sciData[3]
                if not newProgramme:
                    ch -= 1
                #Get first uscount of current package
                scicuruscount = float(struct.unpack('>Q', sciData[4:12])[0]) / specificParam.internalFreq
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
                    scicuruscount = float(struct.unpack('>Q', sciData[27 + 11 * ie:35 + 11 * ie])[0]) / specificParam.internalFreq
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
                sciPattern = re.compile(specificParam.pattrenRef['hexSciNew'], re.I)
                telPattern = re.compile(specificParam.pattrenRef['hexTelNew'], re.I)
            else:
                sciPattern = re.compile(specificParam.pattrenRef['hexSci'], re.I)
                telPattern = re.compile(specificParam.pattrenRef['hexTel'], re.I)
            #Find positions of all science and telemetry data packs
            data = lines[0]
            sciPackPos = basic.findPackPos(data, sciPattern)
            telPackPos = basic.findPackPos(data, telPattern)
            #Get all science and telemetry data packs
            sciPackData = []
            telPackData = []
            for idata in range(len(sciPackPos)):
                sciPackData.append(data[sciPackPos[idata]:sciPackPos[idata] + 512 * 3])
            for idata in range(len(telPackPos)):
                telPackData.append(data[telPackPos[idata]:telPackPos[idata] + 512 * 3])
            if doDeduplicate:
                #Deduplicate data
                sciPackData, telPackData = basic.deduplicate(sciPackData, telPackData)

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
                            if not basic.crcCheck(buf[0:510], bufcrc):
                                #Remove the buffer if CRC check not passed
                                crcError += 1
                                del lineBuffer[lineBuffer.index(buf)]
                            else:
                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                for it in range(7):
                                    utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                    uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                    for ich in range(4):
                                        curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                        tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                        tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                        iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                        bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * specificParam.internalResistance)
                                del lineBuffer[lineBuffer.index(buf)]
                            break
                    #CRC check of current science pack
                    if not basic.crcCheck(lineFloat[:502], lineFloat[502:504]):
                        crcError += 1
                        continue
                else:
                    #CRC check of current science pack
                    if not basic.crcCheck(lineFloat[:510], lineFloat[510:512]):
                        crcError += 1
                        continue
                #Channel of currrent event
                ch = lineFloat[3]
                #Get first uscount of current package
                scicuruscount = float(sum([lineFloat[4 + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq
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
                    scicuruscount = float(sum([lineFloat[27 + 11 * ie + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq
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
                            if not basic.crcCheck(buf[0:510], bufcrc):
                                #Remove the buffer if CRC check not passed
                                crcError += 1
                                del lineBuffer[lineBuffer.index(buf)]
                            else:
                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                for it in range(7):
                                    utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                    uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                    for ich in range(4):
                                        curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                        tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                        tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                        vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                        iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                        bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * specificParam.internalResistance)
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
                    if not basic.crcCheck(lineFloat[:510], temp):
                        crcError += 1
                        continue
                else:
                    #CRC check of current telemetry pack
                    if not basic.crcCheck(lineFloat[:496], lineFloat[496:498]):
                        crcError += 1
                        continue
                for it in range(7):
                    utc.append(float(sum([lineFloat[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                    uscount.append(float(sum([lineFloat[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                    for ich in range(4):
                        curtemp = float(lineFloat[23 + 2 * ich + 70 * it] * 256 + lineFloat[24 + 2 * ich + 70 * it])
                        tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                        curtemp = float(lineFloat[31 + 2 * ich + 70 * it] * 256 + lineFloat[32 + 2 * ich + 70 * it])
                        tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                        vMon[ich].append(float(lineFloat[39 + 2 * ich + 70 * it] * 256 + lineFloat[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                        iMon[ich].append(float(lineFloat[47 + 2 * ich + 70 * it] * 256 + lineFloat[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                        bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * specificParam.internalResistance)

        #--------------------------------------------------Non-hexprint(normal) format files--------------------------------------------------
        else:
            #For some files, ignore CRC check to ensure there is still telemetry data for further processing
            if filenameNoPath in specificParam.ignoreCrcRef:
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
                                            if not basic.crcCheck(bufMatch[0:510], bufcrc):
                                                #Remove the buffer if CRC check not passed
                                                crcError += 1
                                                del lineBuffer[lineBuffer.index(buf)]
                                            else:
                                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                                if isCi == 2:
                                                    for it in range(7):
                                                        utc[nScan].append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount[nScan].append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                                        for ich in range(4):
                                                            curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                                            tempSipm[ich][nScan].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                                            tempAdc[ich][nScan].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            vMon[ich][nScan].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich][nScan].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich][nScan].append(vMon[ich][nScan][-1] - iMon[ich][nScan][-1] * specificParam.internalResistance)
                                                else:
                                                    for it in range(7):
                                                        utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                                        for ich in range(4):
                                                            curtemp = float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it])
                                                            tempSipm[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            curtemp = float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it])
                                                            tempAdc[ich].append((curtemp - 4096) / 16.0 if curtemp > 2048 else curtemp / 16.0)
                                                            vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * specificParam.internalResistance)
                                                del lineBuffer[lineBuffer.index(buf)]
                                            break
                                    #CRC check of current science pack
                                    if not basic.crcCheck(lineFloat[:502], lineFloat[502:504]):
                                        crcError += 1
                                        continue
                                else:
                                    #CRC check of current science pack
                                    if not basic.crcCheck(lineFloat[:510], lineFloat[510:512]):
                                        crcError += 1
                                        continue
                                #Channel of currrent event
                                ch = lineFloat[il + 3]
                                #Get first uscount of current package
                                scicuruscount = float(sum([lineFloat[il + 4 + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq
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
                                    scicuruscount = float(sum([lineFloat[il + 27 + 11 * ie + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq
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
                                            if not basic.crcCheck(bufMatch[0:510], bufcrc):
                                                #Remove the buffer if CRC check not passed
                                                crcError += 1
                                                del lineBuffer[lineBuffer.index(buf)]
                                            else:
                                                #If CRC correct, extract the data of previous data pack and delete matched pack from buffer
                                                if isCi == 2:
                                                    for it in range(7):
                                                        utc[nScan].append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount[nScan].append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                                        for ich in range(4):
                                                            tempSipm[ich][nScan].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] > 2048 \
                                                                else tempSipm[ich][nScan].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it]) / 16.0)
                                                            tempAdc[ich][nScan].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] > 2048 \
                                                                else tempAdc[ich][nScan].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it]) / 16.0)
                                                            vMon[ich][nScan].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich][nScan].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich][nScan].append(vMon[ich][nScan][-1] - iMon[ich][nScan][-1] * specificParam.internalResistance)
                                                else:
                                                    for it in range(7):
                                                        utc.append(float(sum([buf[3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                                        uscount.append(float(sum([buf[15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                                        for ich in range(4):
                                                            tempSipm[ich].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it] > 2048 \
                                                                else tempSipm[ich].append(float(buf[23 + 2 * ich + 70 * it] * 256 + buf[24 + 2 * ich + 70 * it]) / 16.0)
                                                            tempAdc[ich].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                                if buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it] > 2048 \
                                                                else tempAdc[ich].append(float(buf[31 + 2 * ich + 70 * it] * 256 + buf[32 + 2 * ich + 70 * it]) / 16.0)
                                                            vMon[ich].append(float(buf[39 + 2 * ich + 70 * it] * 256 + buf[40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                                            iMon[ich].append(float(buf[47 + 2 * ich + 70 * it] * 256 + buf[48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * specificParam.internalResistance)
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
                                    if not basic.crcCheck(lineFloat[0:510], temp):
                                        crcError += 1
                                        continue
                                else:
                                    #CRC check of current telemetry pack
                                    if not basic.crcCheck(lineFloat[0:496], lineFloat[496:498]) and not ignoreCrc:
                                        crcError += 1
                                        continue
                                for it in range(7):
                                    if isCi == 2:
                                        utc[nScan].append(float(sum([lineFloat[il + 3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                        uscount[nScan].append(float(sum([lineFloat[il + 15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                        for ich in range(4):
                                            tempSipm[ich][nScan].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] > 2048 \
                                                else tempSipm[ich][nScan].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it]) / 16.0)
                                            tempAdc[ich][nScan].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] > 2048 \
                                                else tempAdc[ich][nScan].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it]) / 16.0)
                                            vMon[ich][nScan].append(float(lineFloat[il + 39 + 2 * ich + 70 * it] * 256 + lineFloat[il + 40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                            iMon[ich][nScan].append(float(lineFloat[il + 47 + 2 * ich + 70 * it] * 256 + lineFloat[il + 48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                            bias[ich][nScan].append(vMon[ich][nScan][-1] - iMon[ich][nScan][-1] * specificParam.internalResistance)
                                    else:
                                        utc.append(float(sum([lineFloat[il + 3 + 70 * it + ius] * 256 ** (3 - ius) for ius in range(4)])))
                                        uscount.append(float(sum([lineFloat[il + 15 + 70 * it + ius] * 256 ** (7 - ius) for ius in range(8)])) / specificParam.internalFreq)
                                        for ich in range(4):
                                            tempSipm[ich].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it] > 2048 \
                                                else tempSipm[ich].append(float(lineFloat[il + 23 + 2 * ich + 70 * it] * 256 + lineFloat[il + 24 + 2 * ich + 70 * it]) / 16.0)
                                            tempAdc[ich].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] - 4096) / 16.0) \
                                                if lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it] > 2048 \
                                                else tempAdc[ich].append(float(lineFloat[il + 31 + 2 * ich + 70 * it] * 256 + lineFloat[il + 32 + 2 * ich + 70 * it]) / 16.0)
                                            vMon[ich].append(float(lineFloat[il + 39 + 2 * ich + 70 * it] * 256 + lineFloat[il + 40 + 2 * ich + 70 * it]) / 4096.0 * 3.3 * 11.0)
                                            iMon[ich].append(float(lineFloat[il + 47 + 2 * ich + 70 * it] * 256 + lineFloat[il + 48 + 2 * ich + 70 * it]) / 4096.0 * 3.3)
                                            bias[ich].append(vMon[ich][-1] - iMon[ich][-1] * specificParam.internalResistance)


    #------------------------------------------------------------Readout data preliminary processing------------------------------------------------------------
    if not specificParam.deduplicateFile:
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
            if filenameNoPath in specificParam.cutFileRef:
                timeCut = specificParam.cutFileRef[filenameNoPath]
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
                    curcorr = basic.tempBiasCorrection(tempSipm[:, isc], bias[:, isc])[0]
                    curenergy = None
                    for ich in range(4):
                        curenergy = basic.ecCorrection(np.array(amp[ich][isc]) * curcorr[ich], True, ich)
                        if isinstance(energyCut, list):
                            q1 = (np.array(curenergy) >= energyCut[0]) * (np.array(curenergy) <= energyCut[1])
                        else:
                            q1 = np.array(curenergy) >= energyCut
                        uscountEvt[ich][isc] = np.array(uscountEvt[ich][isc])[q1]
                        amp[ich][isc] = np.array(amp[ich][isc])[q1]
                        if isBin:
                            sciNum[ich][isc] = np.array(sciNum[ich][isc])[q1]
            else:
                curcorr = basic.tempBiasCorrection(tempSipm, bias)[0]
                curenergy = None
                for ich in range(4):
                    curenergy = basic.ecCorrection(np.array(amp[ich]) * curcorr[ich], True, ich)
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
    print('Data readout of \"' + filenameNoPath + '\" complete')

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

#****************************************************************************************************************************************************************
#***************************************************************Basic data process functions**************************************************************
#****************************************************************************************************************************************************************

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
        if not basic.isChannel(channel):
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
                if filename in specificParam.unusedTelRef and curTelNum in specificParam.unusedTelRef[filename]:
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
                if filename in specificParam.unusedSciRef and curSciNum in specificParam.unusedSciRef[filename]:
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
                    curcorr = basic.tempBiasCorrection(tempCorr[:, isci], biasCorr[:, isci])[0]
                    curenergy = None
                    for ich in range(4):
                        curenergy = basic.ecCorrection(np.array(ampFinal[ich][isci]) * curcorr[ich], True, ich)
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
        refVaules = specificParam.biasRef['bias_set']
        bound = specificParam.biasRef['bound']
        if filename in specificParam.unusedRef:
            #Skip unused values in some files
            refVaules = specificParam.unusedRef[filename]
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
    specTime, xdata = np.histogram(timeCorrect, bins = 250, range=[np.min(timeCorrect), np.min([np.max(timeCorrect), specificParam.liveTimeCut])])
    xdata = (xdata[:-1] + xdata[1:]) / 2
    C = specificParam.deadTimeMCU
    center = 0.
    sigma = 0.
    if rateStyle == 's':
        #512-pack dead time correction
        qFit = (specTime > 0) * (xdata > C)
        result = basic.doLogLinFit(xdata[qFit], specTime[qFit], odr = odr, yerror = basic.gehrelsErr(specTime[qFit]))
        a = result['fit_a']
        b = result['fit_b']
        bErr = result['fit_b_err']
        rateAll = b
        rateAllErr = bErr
    else:
        #4096-pack dead time correction
        qFit = (specTime > 0) * (xdata > 43 * C)
        #Gaussian pre-fit for fit range correction
        res = basic.doFitPeak(xdata[qFit], specTime[qFit], odr = False, bkg = False)
        center = res['peak_center']
        sigma = res['peak_sigma']
        #Fit range correction
        specTime, xdata = np.histogram(timeCorrect, bins = 100, range=[np.max([43 * C, center - 6. * sigma]), center + 6. * sigma])
        xdata = (xdata[:-1] + xdata[1:]) / 2
        qFit = (specTime > 0) * (xdata > 43 * C) * (xdata >= center - 4. * sigma) * (xdata <= center + 4. * sigma)
        result = basic.doConvExpFit(xdata[qFit], specTime[qFit], (1., C), odr = odr, yerror = basic.gehrelsErr(specTime[qFit]))
        a = result['fit_a']
        b = result['fit_b']
        bErr = result['fit_b_err']
        rateAll = 1 / b
        rateAllErr = bErr / b ** 2
    if plot:
        if rateStyle == 's':
            fitSpec = np.exp(basic.logLinFunction([a, b], xdata))
            qPlot = xdata > C
        else:
            fitSpec = np.exp(basic.convExpFunction([a, b], xdata))
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

def fileOutput(file, sciData, telData, scanData = {}, isCi = 0, isScan = False, rateStyle = '', newProgramme = False):

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
        with open(parameters.outputConfig['02'], 'r') as fin:
            config = json.load(fin)
    except:
        print('fileOutput: unable to read file output configuration from file \"' + parameters.outputConfig['02'] + '\"')
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
