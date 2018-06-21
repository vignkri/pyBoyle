#!/usr/bin/env python

"""
Analysis Tools

A collection of tools that are useful
for performing analysis on the dataset.
"""

import numpy as np
from scipy import interpolate


def interpolateData(data, time_col, time_step=1):
    """Interpolate the data

    Wrapper to interpolate the result dataset from the
    simulation manager for computing the sparse
    output from the dense output.

    PARAMETERS
    ----------
    data : numpy.array ::
        The numpy array values for the data output from
        the simulation manager.

    time_col : int ::
        The index integer for the column where the time
        data is stored in the array.

    time_step : int ::
        The integer equal to the interpolation frequency
        of the data. Default is set to 1-unit which is
        equal to 1-Hour for the simulation output.

    RETURNS
    -------
    out_ : numpy.array ::
        The interpolated numpy grid data after interpolation
        using the scipy.interpolate.griddata function.
    """
    time = data[:, time_col]
    start_time = time.astype(int).min()
    end_time = time.astype(int).max()
    tarr = np.arange(start_time, end_time, time_step)
    out_ = interpolate.griddata(time, data, tarr, method="linear")
    return out_


def computeBMP(arr, multiplier=None):
    """Compute BMP from the concentration data

    BioMethane Potential of the substrate feeds are computed
    from the debug dataset using the negative yields of the
    substrate. The negative yields of the substrate are
    then matrix multiplied with constants for the different
    components to compute the BMP before correction factors
    for the reactor volume are applied.

    PARAMETERS
    ----------
    arr : numpy.array ::
        Numpy array of the parameters used for the computation
        of the negative yields from the simulation manager.

    multiplier : numpy.array ::
        Numpy array of the multipliers to compute the BMP. The
        default multiplier is set to None to utilise a standard
        constants array.

    RETURNS
    -------
    numpy.array ::
        The numpy array of shape (n, 1) containing the BMP values
        before correcting for the reactor volume.
    """
    if not multiplier:
        multiplier = np.array([0.5*373.5, 373.5, 0.8*572.6, 572.6,
                               1014, 1014, 530, 636.6, 714, 373.5])
    # --
    return (arr) @ (multiplier)
