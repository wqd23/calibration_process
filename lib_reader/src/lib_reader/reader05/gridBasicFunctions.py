# -*- coding:utf-8 -*-

"""
GRID Basic Functions Library
----------

Basic functions for all versions of GRID data processing and fit\n
`v1.0.0` for GRID calibration result analysis, by ghz, ydx and cjr\n
"""

#****************************************************************************************************************************************************************
#**********************************************************************Import packages**********************************************************************
#****************************************************************************************************************************************************************

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
import json
import datetime
from copy import copy

#****************************************************************************************************************************************************************
#**********************************************************Basic data process and fit functions**********************************************************
#****************************************************************************************************************************************************************

#****************************************************************************************************************************************************************
#********************************************************************Auxiliariry functions*******************************************************************
#****************************************************************************************************************************************************************

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

def getSpectrum(amp, nbins = 65536., singlech = False, specRange = [0., 65536.], specEdges = None, binWidth = None, adcMax = 65536.):

    """
    Function for getting the spectrum of the amplitude data

    Parameters
    ----------
    amp : array-like,
        amplitude of all 4 channels\n
    nbins : int, optional
        number of bins, 0 < nbins <= adcMax\n
    singlech : boolean, optional
        indicates whether the spectrum is single-channeled\n
    specRange : list, optional
        list for upper and lower bounds for the spectrum, in the form of `[lower, upper]` with 0 < lower < upper, or \
            `[[lower, upper], [lower, upper], [lower, upper], [lower, upper]]` spectrum ranges are channel-distinct\n
    specEdges : list or array-like or None, optional
        pre-defined spectrum edges, if not None, `specEdges` will be used instead of `nbins` and `specRange`
    binWidth : int or None, optional
        bin width to be used in the spectrum, within [1, adcMax], or None if not used\n
    adcMax : float, optional
        maximum ADC value for current data pack version, which would usually be 'adcMax' in the version-distinct 'gridParametersXX.py' with XX being data \
            pack version

    Returns
    ----------
    spectrum : array-like or list of arrays,
        the corresponding spectrum(s) of the input\n
    x : array-like or list of arrays,
        bin centers of the spectrum(s)\n

    Raises
    ----------
    Exception
        when nbins is not within range `[1, adcMax]`, or specRange is not in correct form\n
    """

    sameRange = True
    if specEdges is None:
        if binWidth is None:
            if nbins <= 0 or nbins > adcMax:
                raise Exception('getSpectrum: parameter \'nbins\' out of range [1-' + str(adcMax) + ']')
        else:
            if binWidth <= 0 or binWidth > adcMax:
                raise Exception('getSpectrum: parameter \'binWidth\' out of range [1-' + str(adcMax) + ']')
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

def checkBinSpecification(filename, binSpecified, wbinSpecified, source = ''):

    """
    Auxiliary function for checking whether the nbins and bin width value of current file (and corresponding source) has been specified in internal parameters \
        (`gridParametersCommon.py`)
    
    Parameters
    ----------
    filename : str, 
        name of current file to check\n
    binSpecified : boolean or None, 
        `binFileSpecified` in `GridDataProcessor.py`, to mark whether all data files have corresponding specified nbins values\n
    wbinSpecified : boolean or None, 
        `wbinFileSpecified` in `GridDataProcessor.py`, to mark whether all data files have corresponding specified bin width values\n
    source : str, optional
        source of corresponding data file\n
    """

    filenameNoPath = extractFilename(filename)
    if filenameNoPath in parameters.nbinsFileRef:
        curbinSpecified = False
        if isinstance(parameters.nbinsFileRef[filenameNoPath], int):
            curbinSpecified = True
        elif source != '' and source in parameters.nbinsFileRef[filenameNoPath]:
            curbinSpecified = True
        binSpecified = curbinSpecified if binSpecified is None else curbinSpecified and binSpecified
    elif filenameNoPath in parameters.binWidthFileRef:
        curwbinSpecified = False
        if isinstance(parameters.binWidthFileRef[filenameNoPath], float):
            curwbinSpecified = True
        elif source != '' and source in parameters.binWidthFileRef[filenameNoPath]:
            curwbinSpecified = True
        wbinSpecified = curwbinSpecified if wbinSpecified is None else curwbinSpecified and wbinSpecified
    return binSpecified, wbinSpecified

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
    return crcCalculated == crcData

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

#****************************************************************************************************************************************************************
#************************************************************************Basic I/O part***********************************************************************
#****************************************************************************************************************************************************************

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

#****************************************************************************************************************************************************************
#***********************************************************************Basic plot part***********************************************************************
#****************************************************************************************************************************************************************

def plotRawData(filename, amp, nbins, corr, time, singlech = False, channel = -1, rateStyle = '', rateAll = None, doCorr = True, utc = None, ecCorr = False, \
    bkg = False, bkgAmp = [], bkgTime = -1.0, bkgRate = None, binWidth = None, ver = '02'):
    
    """
    Function for plotting the processed, unfitted data

    Parameters
    ----------
    filename : str,
        name of raw data file\n
    amp : array of array or array-like,
        ADC amplitude data\n
    nbins : int,
        number of bins to be used in the spectrum, within [1, adcMax]\n
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
    bkg : boolean, optional
        indicates whether the corresponding background is given. The corresponding total measurement times(source and background) should be \
            given(bkgTime) along with the background amplitude\n
    bkgAmp : array-like or array of arrays, optional
        TEMP-BIAS CORRECTED amplitude of background data\n
    bkgTime : float, optional
        total measurement time of background data\n
    bkgRate : float or list, optional
        correct count rate of all background spectrum in all 4 channels calculated when reading data, or count rates for 4 seperate channels, only \
            used when rateStyle is 's' or 'p'\n
    binWidth : int or None, optional
        bin width to be used in the spectrum, within [1, adcMax], or None if not used\n
    ver : str, optional
        data pack version, currently supporting '02' for GRID02 and GRID02-based, and '03' for GRID03 and GRID03-based\n

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
        specRange = [0., parameters.adcMax[ver] * corr[channel]]
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech, channel, ver = ver)
            specRange = [0., parameters.energyPlotUpper[ver]]
        specWidth = specRange[1] - specRange[0]
        spectrumRaw, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
        #Divide bin count by bin width, so that the integral is uniformed to total count
        if binWidth is None:
            spectrumRaw = spectrumRaw * nbins / specWidth
        else:
            spectrumRaw = spectrumRaw / binWidth
    #Multiple(all-4) channel plots
    else:
        dataCorr = []
        specRange = []
        specWidth = []
        for ich in range(4):
            amp[ich] = np.array(amp[ich])
            dataCorr.append(amp[ich] * corr[ich])
            specRange.append([0., parameters.adcMax[ver] * corr[ich]])
            specWidth.append(specRange[ich][1] - specRange[ich][0])
        dataCorr = np.array(dataCorr, dtype = object)
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech, ver = ver)
            specRange = [0., parameters.energyPlotUpper[ver]]
            specWidth = specRange[1] - specRange[0]
        spectrumRaw, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
        for ich in range(4):
            #Divide bin count by bin width, so that the integral is uniformed to total count
            if binWidth is None:
                if isinstance(specWidth, list):
                    spectrumRaw[ich] = spectrumRaw[ich] * nbins / specWidth[ich]
                else:
                    spectrumRaw[ich] = spectrumRaw[ich] * nbins / specWidth
            else:
                spectrumRaw[ich] = spectrumRaw[ich] / binWidth

    #Count rate correction
    rateFactor = 1.
    spectrumCorr = []
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
        if singlech:
            if isinstance(rateAll, list):
                spectrumCorr = spectrumRaw * rateFactor[ich]
            else:
                spectrumCorr = spectrumRaw * rateFactor
        else:
            for ich in range(4):
                if isinstance(rateAll, list):
                    spectrumCorr.append(spectrumRaw[ich] * rateFactor[ich])
                else:
                    spectrumCorr.append(spectrumRaw[ich] * rateFactor)
    else:
        if singlech:
            spectrumCorr = spectrumRaw / time
        else:
            for ich in range(4):
                spectrumCorr.append(spectrumRaw[ich] / time)

    #External background data spectrum and correction
    bkgRateFactor = 1.
    bkgSpectrum = []
    spectrum = []
    if bkg:
        if len(bkgAmp) == 0 or rateStyle == '' and bkgTime < 0:
            raise Exception('fitSpectrum: background data not given')
        if rateStyle == '':
            bkgRateFactor = 1. / bkgTime
        else:
            if isinstance(bkgRate, list):
                bkgRateFactor = []
                for ich in range(4):
                    bkgRateFactor.append(bkgRate[ich] / float(len(bkgAmp[ich])))
            else:
                bkgCountAll = 0.0
                for ich in range(4):
                    bkgCountAll += float(len(bkgAmp[ich]))
                bkgRateFactor = bkgRate / bkgCountAll
        if singlech:
            bkgSpectrum, _ = getSpectrum(bkgAmp[channel], nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
            #Divide bin count by bin width, so that the integral is uniformed to total count
            if binWidth is None:
                bkgSpectrum = bkgSpectrum * nbins / specWidth
            else:
                bkgSpectrum = bkgSpectrum / binWidth
            #Subtraction of external background
            if isinstance(bkgRate, list):
                spectrum = spectrumCorr - bkgSpectrum * bkgRateFactor[channel]
                bkgSpectrum = np.array(bkgSpectrum) * bkgRateFactor[channel]
            else:
                spectrum = spectrumCorr - bkgSpectrum * bkgRateFactor
                bkgSpectrum = np.array(bkgSpectrum) * bkgRateFactor
        else:
            bkgSpectrum, _ = getSpectrum(bkgAmp, nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
            for ich in range(4):
                #Divide bin count by bin width, so that the integral is uniformed to total count
                if binWidth is None:
                    if isinstance(specWidth, list):
                        bkgSpectrum[ich] = bkgSpectrum[ich] * nbins / specWidth[ich]
                    else:
                        bkgSpectrum[ich] = bkgSpectrum[ich] * nbins / specWidth
                else:
                    bkgSpectrum[ich] = bkgSpectrum[ich] / binWidth
            #Subtraction of external background
            if isinstance(bkgRate, list):
                for ich in range(4):
                    spectrum.append(spectrumCorr[ich] - bkgSpectrum[ich] * bkgRateFactor[ich])
                    bkgSpectrum[ich] = np.array(bkgSpectrum[ich]) * bkgRateFactor[ich]
            else:
                for ich in range(4):
                    spectrum.append(spectrumCorr[ich] - bkgSpectrum[ich] * bkgRateFactor)
                    bkgSpectrum[ich] = np.array(bkgSpectrum[ich]) * bkgRateFactor
            spectrum = np.array(spectrum, dtype = object)
    else:
        spectrum = spectrumCorr

    fig = plt.figure(figsize = (12, 8))
    #Single channel plots
    if singlech:
        ich = channel
        gs = gridspec.GridSpec(1, 1, wspace = 0.5, hspace = 0.2, left = 0.13, right = 0.95)
        ax = fig.add_subplot(gs[0])
        plt.step(x, spectrum, where = 'mid', label = 'raw data', zorder = 1)
        if parameters.plotSpecLogScale:
            if ecCorr:
                ax.set_xlim([parameters.logScaleBegin[ver], parameters.energyPlotUpper[ver]])
            else:
                ax.set_xlim([parameters.logScaleBegin[ver], parameters.adcMax[ver]])
            ax.set_xscale('log')
            ax.set_yscale('log')
        else:
            if ecCorr:
                ax.set_xlim([0., parameters.energyPlotUpper[ver]])
            else:
                ax.set_xlim([0., parameters.adcMax[ver]])
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
            plt.step(x[ich], spectrum[ich], where = 'mid', label = 'raw data', zorder = 1)
            if parameters.plotSpecLogScale:
                if ecCorr:
                    ax.set_xlim([parameters.logScaleBegin[ver], parameters.energyPlotUpper[ver]])
                else:
                    ax.set_xlim([parameters.logScaleBegin[ver], parameters.adcMax[ver]])
                ax.set_xscale('log')
                ax.set_yscale('log')
            else:
                if ecCorr:
                    ax.set_xlim([0., parameters.energyPlotUpper[ver]])
                else:
                    ax.set_xlim([0., parameters.adcMax[ver]])
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
                            plt.title('Light curve of file ' + filename + ', run #' + str(irun) + '-' + str(splitcount) + titleAppend)
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
                    plt.title('Light curve of file ' + filename + ', run #' + str(irun) + titleAppend)
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
                plt.title('Light curve of file ' + filename + titleAppend)
            else:
                plt.title('Light curve of file ' + filename + ', scan #' + str(curNum) + titleAppend)
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
                                        ax.set_title('Light curve of file ' + filename + ', run#' + str(irun) + '-' + str(splitcount) + titleAppend)
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
                                    ax.set_title('Light curve of file ' + filename + ', run#' + str(irun) + titleAppend)
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
                    ax.set_title('Light curve of file ' + filename + titleAppend)
                else:
                    ax.set_title('Light curve of file ' + filename + ', scan #' + str(curNum) + titleAppend)
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

#****************************************************************************************************************************************************************
#**************************************************Basic calibration (temperature, bias and EC) part*************************************************
#****************************************************************************************************************************************************************
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

def ecCorrection(amp, singlech = False, channel = -1, ecFiles = parameters.ecConfig, dataCorr = True, ver = '02'):
    #TBD: Add error to corrected energy values
    """
    Function for calculating the EC calibrated values of the spectrum's x-axis(energy)

    Parameters
    ----------
    amp : array of arrays or array-like,
        the input ADC values\n
    singlech : boolean, optional
        indicates whether the input data is single-channeled\n
    channel : int, optional
        the channel number in range [0-3]\n
    ecFiles : dict of list, optional
        lookup table for files with EC correction coeficients, with keys being data pack versions and values being list of file names\n
    dataCorr : boolean, optional
        indicates whether the input is amplitude data, False when input is spectrum fit range to be corrected\n
    ver : str, optional
        data pack version, currently supporting '02' for GRID02 and GRID02-based, and '03' for GRID03 and GRID03-based\n

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
            with open(ecFiles[ver][channel], 'r') as fin:
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

def tempBiasCorrection(temp, bias, doCorr = True, coefFile = parameters.tempBiasConfig, ver = '02', style = 'sci'):

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
    coefFile : dict of str, optional
        lookup table for files with temp-bias correction coeficients, with keys being data pack versions and values being the filenames\n
    ver : str, optional
        data pack version, currently supporting '02' for GRID02 and GRID02-based, and '03' for GRID03 and GRID03-based\n
    style : str, optional
        fit function style for surface fit, with 'sci' being scientific style (physical temperature-bias model) and 'eng' being engineering style (approximated \
            2D-quadratic model)\n

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

    def tempbiasFunctionInternal(param, x, y, xdiff = 0, ydiff = 0, style = 'sci'):

        """
        Auxiliary function to calculate the value of 2D temperature-bias response function and its (1st) derivative, 
        with the formula being\n
        `f(x, y) = G0 * (y - k * x - V0) ^ 2 * (-x ^ 2 + b * x + c)`\n
        for style `sci`, \n
        or 2-d quadratic function and its (1st) derivative, 
        with the formula being\n
        `f(x, y) = p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ^ 2 + p02 * y ^ 2 + p21 * x ^ 2 * y + p12 * x * y ^ 2 + p22 * x ^ 2 * y ^ 2 + p30 + x ^ 3 + \
            p31 * x * 3 * y + p40 * x ^ 4`\n
        for style `eng`

        Parameters
        ----------
        param : list, 
            parameters of the function, in the form of `[G0, k, V0, b, c]`\n
        x : float or array-like,
            value(s) of input x\n
        y : float or array-like,
            value(s) of input y\n
        xdiff : int, optional
            partial derivative order of x, default 0 for `f(x, y)`, and 1 for `df/dx`\n
        ydiff : int, optional
            partial derivative order of y, default 0 for `f(x, y)`, and 1 for `df/dy`\n
        style : str, optional
            fit function style for surface fit, with 'sci' being scientific style (physical temperature-bias model) and 'eng' being engineering style \
                (approximated 2D-quadratic model)\n

        Returns
        ----------
        result : float or array-like,
            values of temperature-bias response function corresponding to input x and given style\n

        Raises
        ----------
        Exception
            when xdiff or ydiff is not 0 or 1, or `style` given is not an available style
        """

        result = None
        if style == 'sci':
            [G0, k, V0, b, c] = param
            if xdiff == 0:
                if ydiff == 0:
                    result = G0 * (y - k * x - V0) ** 2 * (-x ** 2 + b * x + c)
                elif ydiff == 1:
                    result = 2 * G0 * (y - k * x - V0) * (-x ** 2 + b * x + c)
                else:
                    raise Exception('tempbiasFunctionInternal: illegal derivative order ' + str(ydiff))
            elif xdiff == 1:
                if ydiff == 0:
                    result = G0 * (y - k * x - V0) * ((y - k * x - V0) * (-2 * x + b) - k * (-x ** 2  + b * x + c))
                elif ydiff == 1:
                    result = G0 * (2 * (y - k * x - V0) * (-2 * x + b) - k * (-x ** 2  + b * x + c))
                else:
                    raise Exception('tempbiasFunctionInternal: illegal derivative order ' + str(ydiff))
            else:
                raise Exception('tempbiasFunctionInternal: illegal derivative order ' + str(xdiff))
        elif style == 'eng':
            if xdiff == 0:
                if ydiff == 0:
                    result = param[0] + param[1] * x + param[2] * y + param[3] * x * y + param[4] * x ** 2 + param[5] * y ** 2 + param[6] * x ** 2 * y + \
                        param[7] * x * y ** 2 + param[8] * x ** 2 * y ** 2 + param[9] * x ** 3 + param[10] * x ** 3 * y + param[11] * x ** 4
                elif ydiff == 1:
                    result = param[2] + param[3] * x + 2 * param[5] * y + param[6] * x ** 2 + 2 * param[7] * x * y + 2 * param[8] * x ** 2 * y  + param[10] * x ** 3
                else:
                    raise Exception('tempbiasFunctionInternal: illegal derivative order ' + str(ydiff))
            elif xdiff == 1:
                if ydiff == 0:
                    result = param[1] + param[3] * y + 2 * param[4] * x + 2 * param[6] * x * y + param[7] * y ** 2 + 2 * param[8] * x * y ** 2 + 3 * param[9] * x ** 2 + \
                        3 * param[10] * x ** 2 * y + 4 * param[11] * x ** 3
                elif ydiff == 1:
                    result = param[3] + 2 * param[6] * x + 2 * param[7] * y + 4 * param[8] * x * y + 3 * param[10] * x ** 2
                else:
                    raise Exception('tempbiasFunctionInternal: illegal derivative order ' + str(ydiff))
            else:
                raise Exception('tempbiasFunctionInternal: illegal derivative order ' + str(xdiff))
        else:
            raise Exception('tempbiasFunctionInternal: unknown style \"' + style + '\", available styles: \"sci\", \"eng\"')
        return result

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
            with open(coefFile[style][ver], 'r') as fin:
                tempBiasCoef = json.load(fin)
        except:
            raise Exception('tempBiasCorrection: error reading config file \"' + coefFile[ver] + '\"')

        if style == 'sci':
            #Coefficients of 2-D(correlated) temp-bias correction
            #2-D temp-bias response: G0 * (y - k * x - V0) ** 2 * (-x ** 2 + b * x + c)
            G0 = [tempBiasCoef[ich]['G0'] for ich in range(4)]
            G0Err = [tempBiasCoef[ich]['G0_err'] for ich in range(4)]
            k = [tempBiasCoef[ich]['k'] for ich in range(4)]
            kErr = [tempBiasCoef[ich]['k_err'] for ich in range(4)]
            V0 = [tempBiasCoef[ich]['V0'] for ich in range(4)]
            V0Err = [tempBiasCoef[ich]['V0_err'] for ich in range(4)]
            b = [tempBiasCoef[ich]['b'] for ich in range(4)]
            bErr = [tempBiasCoef[ich]['b_err'] for ich in range(4)]
            c = [tempBiasCoef[ich]['c'] for ich in range(4)]
            cErr = [tempBiasCoef[ich]['c_err'] for ich in range(4)]

            param = [G0[ich], k[ich], V0[ich], b[ich], c[ich]]
            paramErr = [G0Err[ich], kErr[ich], V0Err[ich], bErr[ich], cErr[ich]]
            
        elif style == 'eng':
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
        
        #Correction factors and errors
        corrFactor.append(tempbiasFunctionInternal(param, tempStandard, biasStandard, style = style) / tempbiasFunctionInternal(param, \
            tempAvg, biasAvg, style = style))
        tempStdTerm = tempbiasFunctionInternal(param, tempStandard, biasStandard, style = style) * tempbiasFunctionInternal(param, tempAvg, \
            biasAvg, xdiff = 1, style = style) /tempbiasFunctionInternal(param, tempAvg, biasAvg, style = style) ** 2
        biasStdTerm = tempbiasFunctionInternal(param, tempStandard, biasStandard, style = style) * tempbiasFunctionInternal(param, tempAvg, \
            biasAvg, ydiff = 1, style = style) / tempbiasFunctionInternal(param, tempAvg, biasAvg, style = style) ** 2
        pTermStandard = [1., tempStandard, biasStandard, tempStandard * biasStandard, tempStandard ** 2, biasStandard ** 2, \
            tempStandard ** 2 * biasStandard, tempStandard * biasStandard ** 2, tempStandard ** 2 * biasStandard ** 2, tempStandard ** 3, \
                tempStandard ** 3 * biasStandard, tempStandard ** 4]
        pTermAvg = [1., tempAvg, biasAvg, tempAvg * biasAvg, tempAvg ** 2, biasAvg ** 2, tempAvg ** 2 * biasAvg, tempAvg * biasAvg ** 2, \
            tempAvg ** 2 * biasAvg ** 2, tempAvg ** 3, tempAvg ** 3 * biasAvg, tempAvg ** 4]
        corrErrSum = (tempStdTerm * tempStd) ** 2 + (biasStdTerm * biasStd) ** 2
        for it in range(len(param)):
            corrErrSum += ((pTermStandard[it] * tempbiasFunctionInternal(param, tempAvg, biasAvg, style = style) - pTermAvg[it] * \
                tempbiasFunctionInternal(param, tempStandard, biasStandard, style = style)) / tempbiasFunctionInternal(param, tempAvg, biasAvg, \
                style = style) ** 2 * paramErr[it]) ** 2
        corrErr.append(np.sqrt(corrErrSum))

    return corrFactor, corrErr

#****************************************************************************************************************************************************************
#**************************************************************Basic fit functions and fit part************************************************************
#****************************************************************************************************************************************************************

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
    # set bound for param
    # assume center is in the xdata range, xdata range bigger than 3*sigma(half the peak)
    # param['peak_amplitude'].min, param['peak_amplitude'].max = 0, 10*np.max(ydata)
    param['peak_center'].min, param['peak_center'].max = xdata[0], xdata[-1]
    param['peak_sigma'].min, param['peak_sigma'].max = 0, (xdata[-1]-xdata[0])/3.
    # make sure bkg is positive
    # if bkgForm == 'lin':
    #     param.add('left_point', min=0,vary=True)
    #     param.add('right_point', min=0,vary=True)
    #     param['bk_slope'].set(expr=f"(right_point-left_point)/({xdata[-1]-xdata[0]})")
    #     param['bk_intercept'].set(expr=f"left_point-bk_slope*{xdata[0]}")

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
    `f(x, y) = p00 + p10 * x + p01 * y + p11 * x * y + p20 * x ^ 2 + p02 * y ^ 2 + p21 * x ^ 2 * y + p12 * x * y ^ 2 + p22 * x ^ 2 * y ^ 2 + p30 * x ^ 3 + \
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

def tempBias2DFunction(input, G0, k, V0, b, c):

    """
    Auxiliary function to calculate the value of 2D temperature-bias response function for correlative temperature-bias fit\n
    in the form of\n
    `f(x, y) = G0 * (y - k * x - V0) ^ 2 * (-x ^ 2 + b * x + c) + C0`

    Parameters
    ----------
    input : array-like,
        array containing both x and y data, being `x = input[:, 0]` and `y = input[:, 1]`\n
    param : list, 
        fit parameters in the form of\n
        param = `[G0, k, V0, b, c]`\n
        where the parameters stand for: \n
            'G0':   `G0` term, as the constant defining overall values
            'k':     `k` term, representing the coefficient for breakdown voltage's variation with temperature
            'V0':   `V0` term, representing the reference breakdown voltage value
            'b':      `b` (normalized linear) term, in the Taylor expansion of GAGG light yield's variation with temperature
            'c':      `c` (normalized constant) term, in the Taylor expansion of GAGG light yield's variation with temperature

    Returns
    ----------
    float or array-like,
        values of temperature-bias response function corresponding to input x and y\n
    """

    xdata = input[:, 0]
    ydata = input[:, 1]

    Vov = ydata - k * xdata - V0
    return G0 * Vov ** 2 * (-xdata ** 2 + b * xdata + c)

def residualTempbias2D(param, xdata, ydata, zdata):

    """
    Auxiliary function to calculate the residual of 2D temperature-bias response model for correlative temperature-bias fit\n
    in the form of\n
    `res = zdata - f(xdata, ydata)`\n
    where\n
    `f(x, y) = G0 * (y - k * x - V0) ^ 2 * (-x ^ 2 + b * x + c)`

    Parameters
    ----------
    param : list, 
        fit parameters  in the form of in the form of `[G0, k, V0, b, c]`\n
    xdata : array-like,
        value(s) of input x\n
    ydata : array-like,
        value(s) of input y\n
    zdata : float or array-like,
        value(s) of input z\n

    Returns
    ----------
    float or array-like,
        calculated residual
    """

    #Parameters
    [G0, k, V0, b, c] = param

    return zdata - tempBias2DFunction(np.dstack((xdata, ydata))[0], G0, k, V0, b, c)

def doFitTempbias2D(xdata, ydata, zdata, odr = False, xerror = [], yerror = [], zerror = [], refVal = [], refStd = []):

    """
    Function for fitting 2D temperature-bias response data, namely for temperature-bias response experiments

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
    refVal : list, optional
        reference values for the paramerers, in the form of `[G0, k, V0, b, c]`\n
    refStd : list, optional
        reference standard deviation values for the paramerers, in the form of `[G0, k, V0, b, c]`\n

    Returns
    ----------
    fitResult : dict,
        fit result in the form of a dictionary as\n
        fitResult = {
            'G0':                   `G0` term, as the constant defining overall values,
            'k':                     `k` term, representing the coefficient for breakdown voltage's variation with temperature,
            'V0':                   `V0` term, representing the reference breakdown voltage value,
            'b':                     `b` (normalized linear) term, in the Taylor expansion of GAGG light yield's variation with temperature, 
            'c':                     `c` (normalized constant) term, in the Taylor expansion of GAGG light yield's variation with temperature, 
            'G0_err':             error of `G0` term, as the constant defining overall values,
            'k_err':               error of `k` term, representing the coefficient for breakdown voltage's variation with temperature,
            'V0_err':             error of `V0` term, representing the reference breakdown voltage value,
            'b_err':               error of `b` (normalized linear) term, in the Taylor expansion of GAGG light yield's variation with temperature, 
            'c_err':               error of `c` (normalized constant) term, in the Taylor expansion of GAGG light yield's variation with temperature, 
            'chisquare':         chi-square value to determine goodness of the fit, 
        }
    with fit function as\n
    `f(x, y) = G0 * (y - k * x - V0) ^ 2 * (-x ^ 2 + b * x + c)`
    """
    
    def tempbias2DFunctionInternal(param, x):

        """
        Auxiliary function to calculate the value of 2D temperature-bias response function for correlative temperature-bias fit, used for ODR fit\n
        in the form of\n
        `f(x, y) = G0 * (y - k * x - V0) ^ 2 * (-x ^ 2 + b * x + c)`

        Parameters
        ----------
        param : list, 
            parameters of the function, in the form of `[G0, k, V0, b, c]`\n
        x : array-like,
            array containing both x and y data, being `xdata = x[0]` and `ydata = x[1]`\n

        Returns
        ----------
        float or array-like,
            values of temperature-bias response function corresponding to input x and y\n
        """

        xdata = x[0]
        ydata = x[1]
        [G0, k, V0, b, c] = param
        Vov = ydata - k * xdata - V0
        return G0 * Vov ** 2 * (-xdata ** 2 + b * xdata + c)

    if len(refVal) == 5:
        [G0init, k0, Vbd0, b0, c0] = refVal
    else:
        k0 = 18.9e-3
        Vbd0 = 24.6
        b0 = 0.
        c0 = 1e4
        Vov = parameters.biasStandard - k0 * parameters.tempStandard - Vbd0
        zdataStandard = zdata[np.where((xdata >= parameters.tempStandard - 5 * parameters.tempErrSys) * (xdata <= parameters.tempStandard + \
            5 * parameters.tempErrSys) * (ydata >= parameters.biasStandard - parameters.biasSetStd) * (ydata <= parameters.biasStandard + \
            parameters.biasSetStd))[0][0]]
        G0init = zdataStandard / Vov ** 2 / quadFunction([-1., b0, c0], parameters.tempStandard)
    paramInit = [G0init * 10., k0, Vbd0, b0, c0]
    if len(refStd) == 5:
        [G0Std, k0Std, V0Std, b0Std, c0Std] = refStd
        lowerBound = [max(0, G0init - 6 * G0Std), max(0, k0 - 6 * k0Std), max(0, Vbd0 - 6 * V0Std), b0 - 6 * b0Std, max(0, c0 - 6 * c0Std)]
        upperBound = [G0init + 6 * G0Std, k0 + 6 * k0Std, Vbd0 + 6 * V0Std, b0 + 6 * b0Std, c0 + 6 * c0Std]
    else:
        lowerBound = [1e-5, 1e-3, 23., -500., 1e3]
        upperBound = [100., 50e-3, 25., 500., 1e8]

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
        model = Model(tempbias2DFunctionInternal)
        odrFit = ODR(data, model, paramInit)
        odrFit.set_job(fit_type = 0)
        result = odrFit.run()
        chisq = result.sum_square / (len(zdata) - len(result.beta))
        fitResult = {
                        'G0':             result.beta[0], 
                        'k':               result.beta[1], 
                        'V0':             result.beta[2], 
                        'b':               result.beta[3], 
                        'c':               result.beta[4], 
                        'G0_err':             result.sd_beta[0], 
                        'k_err':               result.sd_beta[1], 
                        'V0_err':             result.sd_beta[2], 
                        'b_err':               result.sd_beta[3], 
                        'c_err':               result.sd_beta[4], 
                        'chisquare':          chisq, 
            }
    else:
        inputData = np.array(np.dstack((xdata, ydata))[0])
        popt, pcov = curve_fit(tempBias2DFunction, inputData, zdata, sigma = zerror, absolute_sigma = True, method = 'trf', p0 = paramInit, \
            bounds = (lowerBound, upperBound))
        param = popt.tolist()
        perr = np.sqrt(np.diag(pcov))
        chisq = sum((residualTempbias2D(param, xdata, ydata, zdata) / zerror) ** 2) / (len(zdata) - len(param))
        fitResult = {
                        'G0':             param[0], 
                        'k':               param[1], 
                        'V0':             param[2], 
                        'b':               param[3], 
                        'c':               param[4], 
                        'G0_err':             perr[0], 
                        'k_err':               perr[1], 
                        'V0_err':             perr[2], 
                        'b_err':               perr[3], 
                        'c_err':               perr[4], 
                        'chisquare':          chisq, 
            }
    paramList = ['G0', 'k', 'V0', 'b', 'c']
    for ip in range(len(paramList)):
        if np.abs((fitResult[paramList[ip]] - lowerBound[ip]) / lowerBound[ip]) < 1e-3:
            print('doFitTempbias2D: param \"' + paramList[ip] + '\" converged at lower bound')
        elif np.abs((fitResult[paramList[ip]] - upperBound[ip]) / upperBound[ip]) < 1e-3:
            print('doFitTempbias2D: param \"' + paramList[ip] + '\" converged at upper bound')
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

#------------------------------------------------------------Spectrum fit function------------------------------------------------------------

#An idea to utilize error of correction factor: Use error propagation:
#   Let C be temp-bias correction factor, f(x) be original distribution, then
#   the corrected distribution is f1(x, C) = f(x / C)
#   Therofore the error propagation can be calculated as
#   sigma_f1 ** 2 = sigma_f ** 2 + (f'(x / C) * x / C ** 2) ** 2 * sigma_C ** 2
#Forget it. Won't work with the amplitude-calibration method currently used.
#Just do Monte-Carlo simulation to describe the error brought by temp-bias corretcion
def peak_fit_bound(xdata, ydata, odr = False, xerror = [], yerror = [], bkg = True, bkgForm = 'lin'):
    import util_lib as util
    return util.peak_fit(xdata, ydata, yerror, bkg_form="lin")
def fitSpectrum(filename, amp, nbins, source, corr, time, fileOutput = False, singlech = False, bkg = False, odr = False, xRange = [],\
    channel = -1, corrErr = [], bkgAmp = [], bkgTime = -1.0, maxiter = 1, bound = 3.0, plot = True, rateStyle = '', rateAll = None, rateAllErr = None,\
    bkgRate = None, bkgRateErr = None, bkgForm = '', poissonErr = False, doCorr = True, lowerLim = -1.0, vLine = [], ecCorr = False, binWidth = None, \
    ver = '02'):

    """
    Function for fitting the spectrum

    Parameters
    ----------
    filename : str,
        name of input data file\n
    amp : array-like or array of arrays,
        ADC amplitude data\n
    nbins : int,
        number of bins to be used in the spectrum, within [1, adcMax]\n
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
        bin width to be used in the spectrum, within [1, adcMax], or None if not used\n
    ver : str, optional
        data pack version, currently supporting '02' for GRID02 and GRID02-based, and '03' for GRID03 and GRID03-based\n

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

    #Read output config file
    config = {}
    try:
        with open(parameters.outputConfig[ver], 'r') as fin:
            config = json.load(fin)
    except:
        raise Exception('fitSpectrum: unable to read file output configuration from file \"' + parameters.outputConfig[ver] + '\"')

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
    if not source in fitRange[ver]:
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
            rangeLim = fitRange[ver][source][channel]
            rangeLim[0] /= rangeCorr[channel]
            rangeLim[1] /= rangeCorr[channel]
            if len(xRange) > 0 and len(xRange[channel]) == 2:
                rangeLim = xRange[channel]
        else:
            rangeLim = fitRange[ver][source]
            for ich in range(4):
                rangeLim[ich][0] /= rangeCorr[ich]
                rangeLim[ich][1] /= rangeCorr[ich]
            if len(xRange) > 0:
                for ich in range(4):
                    if len(xRange[ich]) == 2:
                        rangeLim[ich] = xRange[ich]

    #EC correction for fit range
    if ecCorr:
        rangeLim = ecCorrection(rangeLim, singlech, channel, dataCorr = False, ver = ver)

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
        specRange = [0., parameters.adcMax[ver] * corr[channel]]
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech, channel, ver = ver)
            specRange = [0., parameters.energyPlotUpper[ver]]
        specWidth = specRange[1] - specRange[0]
        spectrumRaw, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
        spectrumStatErr = gehrelsErr(spectrumRaw)
    else:
        specRange = []
        specWidth = []
        dataCorr = []
        for ich in range(4):
            amp[ich] = np.array(amp[ich])
            dataCorr.append(amp[ich] * corr[ich])
            specRange.append([0., parameters.adcMax[ver] * corr[ich]])
            specWidth.append(specRange[ich][1] - specRange[ich][0])
        dataCorr = np.array(dataCorr, dtype = object)
        if ecCorr:
            dataCorr = ecCorrection(dataCorr, singlech, ver = ver)
            specRange = [0., parameters.energyPlotUpper[ver]]
            specWidth = specRange[1] - specRange[0]
        spectrumRaw, x = getSpectrum(dataCorr, nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
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
            bkgSpectrum, _ = getSpectrum(bkgAmp[channel], nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
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
            bkgSpectrum, _ = getSpectrum(bkgAmp, nbins, singlech, specRange, binWidth = binWidth, adcMax = parameters.adcMax[ver])
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
        # result = peak_fit_bound(x[q1], spectrum[q1], odr, yerror = spectrumErr[q1], bkg = (bkgForm != ''), bkgForm = bkgForm)
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
                plotLowerLim = parameters.logScaleBegin[ver]
            lower = center - 6 * sigma
            if lower < plotLowerLim:
                lower = plotLowerLim
            upper = center + 6 * sigma
            if upper > parameters.adcMax[ver] * corr[ich]:
                upper = parameters.adcMax[ver] * corr[ich]
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
                    plotLowerLim = parameters.logScaleBegin[ver]
                lower = center - 6 * sigma
                if lower < plotLowerLim:
                    lower = plotLowerLim
                upper = center + 6 * sigma
                if upper > parameters.adcMax[ver] * corr[ich]:
                    upper = parameters.adcMax[ver] * corr[ich]
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
