# -*- coding:utf-8 -*-
"""Detector resolution model shared by the pipeline and plotting layers."""
import numpy as np


def resolutionFunction(x, a, b, c):
    y = np.sqrt(a * x * x + b * x + c) / x
    y[a * x * x + b * x + c < 0] = 0.0
    return y
