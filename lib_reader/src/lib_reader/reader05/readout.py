# -*- coding:utf-8 -*-
"""Self-contained decoders for the GRID 03B / 05B byte stream.

Extracted from the legacy ``gridBasicFunctions`` / ``gridProcessFunctions03``
modules so the unified reader depends only on these small functions: UDP
science-payload extraction, HK / timeline record decoding (both byte orders) and
the UTC-timestamp fit.  Function bodies are kept byte-for-byte identical to the
legacy implementation.
"""
import re
import struct

import numpy as np
import lmfit
from scipy.odr import ODR, Model, RealData

INTERNAL_FREQ = 100.0e6
UDP_PACK_LEN = 16384
HK_DATA_LEN = 95
TIMELINE_DATA_LEN = 16
MAX_UDP_READOUT = 10000
SPLIT_RUN_TIME = 30.0

PATTERNS = {
    "udp": rb'\x1a\xcf\xfc\x1d.{2}\x11\x19.{2}\x88.{16384}.{1}\x2e\xe9\xc8\xfd',
    "HK_new": rb'\x1a\xcf\xfc\x1d\x00\x00\x00\x00.{4}\x00\x00\x00\x01',
    "HK_03b": b'\\x1a\\xcf\\xfc\\x1d.{2}\\x11\\x19.{2}\\x89.{1}.{52}.{1}\\x2e\\xe9\\xc8\\xfd',
    "timeline": rb'\x1a\xcf\xfc\x1d.{2}\x11\x19.{2}\x90.{16}.{1}\x2e\xe9\xc8\xfd',
    "time_new": rb'\x1a\xcf\xfc\x1d\x00\x00\x00\x00.{4}\x00\x00\x00\x07',
}

# Per-file science time cut (seconds), keyed by file name (legacy cutFileRef).
CUT_FILE_REF = {
    '13.4_ch2_100s_rundata2020-09-05-20-24-52.dat': 12020,
    '15.6_ch0_60s_rundata2020-09-05-20-12-16.dat': 11260,
    '20.9_ch0_30s_rundata2020-09-05-19-50-11.dat': 9900,
    '24.2_ch2_30s_rundata2020-09-05-19-40-52.dat': 9372,
    '34.5_ch3_30s_rundata2020-09-05-19-13-14.dat': 7724,
    '38.3_ch1_40s_rundata2020-09-05-19-03-25.dat': 7120,
    '42.5_ch1_60s_rundata2020-09-05-18-49-11.dat': 6280,
    '49.2_ch3_30s_rundata2020-09-16-13-49-00.dat': 4541,
    '51.4_ch0_30s_rundata2020-09-16-14-18-21.dat': 6300,
    '92.0_ch0_60s_rundata2020-09-17-13-46-38.dat': 12200,
    '92.0_ch2_60s_rundata2020-09-17-13-48-41.dat': 12674,
    '92.0_ch3_60s_rundata2020-09-17-13-51-08.dat': 12571,
    '101.8_ch2_60s_rundata2020-09-17-14-05-56.dat': 13595,
    '113.7_ch0_60s_rundata2020-09-17-14-21-03.dat': 14400,
    '10C_28.0V_4m_rundata2020-09-11-18-46-06.dat': 2500,
    '15C_29.0V_4m_rundata2020-09-11-21-41-17.dat': 6940,
    '20C_26.5V_4m_rundata2020-09-11-23-22-50.dat': 13000,
    '25C_26.5V_4m_rundata2020-09-12-02-06-13.dat': 22800,
    '25C_28.3V_4m_rundata2020-09-12-02-46-06.dat': 25200,
    '30C_27.5V_4m_rundata2020-09-12-05-08-21.dat': 5050,
    '35C_27.5V_4m_rundata2020-09-12-08-09-29.dat': 3030,
    '40C_28.5V_4m_rundata2020-09-12-10-55-33.dat': 3570,
    'jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat': 5025,
    'jly_24p2_ch0_30s_rundata2021-04-29-15-10-52.dat': 4380,
    'jly_31p2_ch3_30s_rundata2021-04-29-14-57-31.dat': 3580,
    'jly_38p2_ch0_30s_rundata2021-04-29-14-38-47.dat': 2460,
    'jly_45p9_ch0_30s_rundata2021-04-29-17-41-58.dat': 6050,
    'jly_45p9_ch3_30s_rundata2021-04-29-17-38-24.dat': 5840,
    'jly_77p0_ch3_30s_rundata2021-04-29-19-00-01.dat': 10740,
}


def findPackPos(data, pattern, getLen=False):
    pos = []
    packLen = []
    for ip in pattern.finditer(data):
        pos.append(ip.start())
        packLen.append(ip.end() - ip.start())
    if getLen:
        return np.array(pos), np.array(packLen)
    else:
        return np.array(pos)


def extractSciRawData(rawData, udpPackPos, maxUdpReadout=-1, lastUdpPos=0, udpPackLen=UDP_PACK_LEN):
    udpPackagesID = []
    sciRawDataList = []
    if maxUdpReadout <= 0:
        posEnd = len(udpPackPos)
    else:
        posEnd = min(lastUdpPos + maxUdpReadout, len(udpPackPos))
    for ipos in range(lastUdpPos, posEnd):
        sciRawDataList.extend(rawData[udpPackPos[ipos] + 11:udpPackPos[ipos] + 11 + udpPackLen])
        udpPackagesID.append(rawData[udpPackPos[ipos] + 8] * 256 + rawData[udpPackPos[ipos] + 9])
        if (len(udpPackagesID) > 2 and (udpPackagesID[-1] - udpPackagesID[-2]) > 100):
            sciRawDataList.clear()
            udpPackagesID.clear()
    udpPackageGap = np.uint16(udpPackagesID[1:]) - np.uint16(udpPackagesID[:-1])
    udpPackageLoss = sum(udpPackageGap) - len(udpPackageGap)
    if udpPackageLoss < 0:
        udpPackageLoss = 0
    sciRawData = bytes(sciRawDataList)
    sciRawDataList.clear()

    return sciRawData


def extractHKData_03b(rawData, hkPackLen=52):
    udpPattern = re.compile(PATTERNS['HK_03b'], re.S)
    udpPos = findPackPos(rawData, udpPattern)

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

    channelLookup = [0, 3, 2, 1]
    for i in np.arange(len(udpPos)):
        HKdata = rawData[udpPos[i] + 12:udpPos[i] + 12 + hkPackLen]
        for ich in np.arange(4):
            bias[channelLookup[ich]].append(struct.unpack('>H', HKdata[2 * ich + 0:2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('>H', HKdata[8 + 2 * ich:8 + 2 * ich + 2])[0])
            curtemp = struct.unpack('>H', HKdata[16 + 2 * ich:16 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('>H', HKdata[32 + 2 * ich:32 + 2 * ich + 2])[0])
        timestamp.append(struct.unpack('>Q', HKdata[24:24 + 8])[0] / INTERNAL_FREQ)

    hkData = {
        'iMon': np.array(iMon) / 2 ** 12 * 2.5 / (1 + 49.9 / 499) / 499 * 1E6,
        'bias': np.array(bias) / 2 ** 12 * 2.5 / (51.1 / (1000 + 51.1)),
        'temp': np.array(temp) / 2 ** 4 * 0.0625,
        'timestamp': np.array(timestamp) * 100.,
        'iSys': np.array(iSys) / 2 ** 12 * 2.5 / (0.05 * 4.7E3 / 100),
    }
    hkData['bias'] = hkData['bias'] - hkData['iMon'] * 499 * 1E-6
    return hkData


def extractHKData_normal(rawData, hkPackLen=HK_DATA_LEN):
    udpPattern = re.compile(PATTERNS['HK_new'], re.S)
    udpPos = findPackPos(rawData, udpPattern)

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

    channelLookup = [0, 1, 2, 3]
    for i in np.arange(len(udpPos)):
        HKdata = rawData[udpPos[i] + 36:udpPos[i] + 36 + hkPackLen]
        for ich in np.arange(4):
            bias[channelLookup[ich]].append(struct.unpack('>H', HKdata[23 + 2 * ich:23 + 2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('>H', HKdata[31 + 2 * ich:31 + 2 * ich + 2])[0])
            curtemp = struct.unpack('>H', HKdata[47 + 2 * ich:47 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('>H', HKdata[63 + 2 * ich:63 + 2 * ich + 2])[0])
        timestamp.append(struct.unpack('>Q', HKdata[55:55 + 8])[0] / INTERNAL_FREQ)

    hkData = {
        'iMon': np.array(iMon) / 2 ** 12 * 2.5 / (1 + 49.9 / 499) / 499 * 1E6,
        'bias': np.array(bias) / 2 ** 12 * 2.5 / (51.1 / (1000 + 51.1)),
        'temp': np.array(temp) / 2 ** 4 * 0.0625,
        'timestamp': np.array(timestamp) * 100.,
        'iSys': np.array(iSys) / 2 ** 12 * 2.5 / (0.05 * 4.7E3 / 100),
    }
    hkData['bias'] = hkData['bias'] - hkData['iMon'] * 499 * 1E-6
    return hkData


def extractHKData(rawData, hkPackLen=HK_DATA_LEN):
    udpPattern = re.compile(PATTERNS['HK_new'], re.S)
    udpPos = findPackPos(rawData, udpPattern)

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

    channelLookup = [0, 1, 2, 3]
    for i in np.arange(len(udpPos)):
        HKdata = rawData[udpPos[i] + 36:udpPos[i] + 36 + hkPackLen]
        for ich in np.arange(4):
            bias[channelLookup[ich]].append(struct.unpack('<H', HKdata[23 + 2 * ich:23 + 2 * ich + 2])[0])
            iMon[channelLookup[ich]].append(struct.unpack('<H', HKdata[31 + 2 * ich:31 + 2 * ich + 2])[0])
            curtemp = struct.unpack('<H', HKdata[47 + 2 * ich:47 + 2 * ich + 2])[0]
            temp[channelLookup[ich]].append(curtemp - 65536 if curtemp > 32768 else curtemp)
            iSys[ich].append(struct.unpack('<H', HKdata[63 + 2 * ich:63 + 2 * ich + 2])[0])

        timestamp.append(struct.unpack('<Q', HKdata[55:55 + 8])[0] / INTERNAL_FREQ)

    hkData = {
        'iMon': np.array(iMon) / 2 ** 12 * 2.5 / (1 + 49.9 / 499) / 499 * 1E6,
        'bias': np.array(bias) / 2 ** 12 * 2.5 / (51.1 / (1000 + 51.1)),
        'temp': np.array(temp) / 2 ** 4 * 0.0625,
        'timestamp': np.array(timestamp) * 100.,
        'iSys': np.array(iSys) / 2 ** 12 * 2.5 / (0.05 * 4.7E3 / 100),
    }
    hkData['bias'] = hkData['bias'] - hkData['iMon'] * 499 * 1E-6
    return hkData


def extractTimelineData(rawData, tlPackLen=TIMELINE_DATA_LEN):
    udpPattern = re.compile(PATTERNS['time_new'], re.S)
    udpPos = findPackPos(rawData, udpPattern)

    tlData = {}
    utc = []
    pps = []
    timestamp = []

    for i in np.arange(len(udpPos)):
        TLdata = rawData[udpPos[i] + 36:udpPos[i] + 36 + tlPackLen]
        utc.append(struct.unpack('>L', TLdata[:4])[0])
        pps.append(struct.unpack('>L', TLdata[4:8])[0])
        timestamp.append(struct.unpack('>Q', TLdata[8:16])[0] / INTERNAL_FREQ)

    tlData = {
        'utc': np.array(utc),
        'pps': np.array(pps),
        'timestamp': np.array(timestamp) * 100.,
    }

    return tlData


def extractTimelineData_03b(rawData, tlPackLen=TIMELINE_DATA_LEN):
    udpPattern = re.compile(PATTERNS['timeline'], re.S)
    udpPos = findPackPos(rawData, udpPattern)

    tlData = {}
    utc = []
    pps = []
    timestamp = []

    for i in np.arange(len(udpPos)):
        TLdata = rawData[udpPos[i] + 11:udpPos[i] + 11 + tlPackLen]
        utc.append(struct.unpack('>L', TLdata[:4])[0])
        pps.append(struct.unpack('>L', TLdata[4:8])[0])
        timestamp.append(struct.unpack('>Q', TLdata[8:16])[0] / INTERNAL_FREQ)

    tlData = {
        'utc': np.array(utc),
        'pps': np.array(pps),
        'timestamp': np.array(timestamp),
    }

    return tlData


def linearFunction(param, x):
    return param[0] * x + param[1]


def doFitLin(xdata, ydata, odr=False, xerror=[], yerror=[]):
    lModel = lmfit.models.LinearModel(prefix='fit_')
    param = lModel.guess(ydata, x=xdata)
    if odr:
        if len(xerror) == 0:
            if len(yerror) == 0:
                data = RealData(xdata, ydata)
            else:
                data = RealData(xdata, ydata, sy=yerror)
        else:
            if len(yerror) == 0:
                data = RealData(xdata, ydata, sx=xerror, fix=np.ones(len(xerror)))
            else:
                data = RealData(xdata, ydata, sx=xerror, sy=yerror, fix=np.ones(len(xerror)))
        model = Model(linearFunction)
        odrFit = ODR(data, model, [param.valuesdict()['fit_slope'], param.valuesdict()['fit_intercept']])
        odrFit.set_job(fit_type=0)
        result = odrFit.run()
        fitResult = {
            'fit_a': result.beta[0],
            'fit_b': result.beta[1],
            'fit_a_err': result.sd_beta[0],
            'fit_b_err': result.sd_beta[1],
        }
    else:
        if len(yerror) == 0:
            result = lModel.fit(ydata, param, x=xdata)
        else:
            result = lModel.fit(ydata, param, x=xdata, weights=1. / np.array(yerror))
        fitResult = {
            'fit_a': result.best_values['fit_slope'],
            'fit_b': result.best_values['fit_intercept'],
            'fit_a_err': result.params['fit_slope'].stderr,
            'fit_b_err': result.params['fit_intercept'].stderr,
        }
    return fitResult


def getUTC(filename, timestampTL, utc, timestamp, plot=False):
    result = doFitLin(timestampTL, utc)
    utcFinal = result['fit_a'] * np.array(timestamp) + result['fit_b']
    utcFinal = np.array(utcFinal, dtype=int)
    return utcFinal
