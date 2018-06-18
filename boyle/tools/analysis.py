#!/usr/bin/env python

"""
Analysis Tools

A collection of tools that are useful
for performing analysis on the dataset.
"""

import numpy as np
from scipy import interpolate


def interpolateData(data, time_col, time_step=1):
    """Interpolate the data"""
    time = data[:, time_col]
    start_time = time.astype(int).min()
    end_time = time.astype(int).max()
    tarr = np.arange(start_time, end_time, time_step)
    out_ = interpolate.griddata(time, data, tarr, method="linear")
    return out_


def computeBMP(arr, multiplier=None):
    """Compute BMP from the concentration data"""
    if not multiplier:
        multiplier = np.array([0.5*373.5, 373.5, 0.8*572.6, 572.6,
                               1014, 1014, 530, 636.6, 714, 373.5])
    # --
    return (arr) @ (multiplier)
