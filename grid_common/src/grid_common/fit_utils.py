# -*- coding:utf-8 -*-
"""Shared spectrum / response helpers.

Extracted verbatim from the legacy ``reader05.gridBasicFunctions`` (via
``lib_reader.reader05.fit_utils``) so the reader package no longer owns them.
Function bodies are unchanged.
"""
import numpy as np


def getSpectrum(amp, nbins=65536., singlech=False, specRange=[0., 65536.], specEdges=None, binWidth=None, adcMax=65536.):
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
                    specch, xch = np.histogram(amp[ich], bins=nbins, range=curRange)
                else:
                    binEdges = np.arange(curRange[0], curRange[1] + binWidth, binWidth)
                    specch, xch = np.histogram(amp[ich], bins=binEdges)
            else:
                specch, xch = np.histogram(amp[ich], bins=specEdges)
            spectrum.append(specch)
            x.append((xch[:-1] + xch[1:]) / 2)
    else:
        if specEdges is None:
            if binWidth is None:
                spectrum, x = np.histogram(amp, bins=nbins, range=(specRange[0], specRange[1]))
            else:
                binEdges = np.arange(specRange[0], specRange[1] + binWidth, binWidth)
                spectrum, x = np.histogram(amp, bins=binEdges)
        else:
            spectrum, x = np.histogram(amp, bins=specEdges)
        x = (x[:-1] + x[1:]) / 2
    return spectrum, x


def gehrelsErr(ydata):
    if np.size(ydata) == 1:
        yerr = np.sqrt(ydata) if ydata >= 5.0 else 1.0 + np.sqrt(ydata + 0.75)
    elif np.size(ydata) > 1:
        yerr = np.sqrt(ydata)
        q = np.where(ydata < 5.0)
        yerr[q] = 1.0 + np.sqrt(ydata[q] + 0.75)
    else:
        yerr = []
    return yerr


def tempBias2DFunction(input, G0, k, V0, b, c):
    xdata = input[:, 0]
    ydata = input[:, 1]

    Vov = ydata - k * xdata - V0
    return G0 * Vov ** 2 * (-xdata ** 2 + b * xdata + c)


def residualTempbias2D(param, xdata, ydata, zdata):
    [G0, k, V0, b, c] = param

    return zdata - tempBias2DFunction(np.dstack((xdata, ydata))[0], G0, k, V0, b, c)


def gaussianFunction(param, x):
    return param[0] * np.exp(-(x - param[1]) ** 2 / (2 * (param[2] ** 2))) / (param[2] * np.sqrt(2 * np.pi))


def quadFunction(param, x):
    return param[0] * x * x + param[1] * x + param[2]
