#!/usr/bin/env python

"""
Analysis Tools

A collection of tools that are useful
for performing analysis on the dataset.
"""

import numpy as np
from scipy import interpolate


def interpolateData(data, time_step=1):
    """Interpolate the data"""
    time = data[:, 1]
    start_time = time.astype(int).min()
    end_time = time.astype(int).max()
    tarr = np.arange(start_time, end_time, time_step)
    out_ = interpolate.griddata(time, data, tarr, method="linear")
    return out_
