#!/usr/bin/env/python
"""
Utility

A collection of utility functions that can
be used by multiple scripts if needed. This is just for
allowing different tools accesses to common functions.
"""

import numpy as np


def load_constants(path):
    """Utility function to load constants files

    Utility function for loading constants files stored
    as plain-text. This file uses '%' as the comment
    delimiter with space delimits for the values in the
    dataset.

    PARAMETERS
    ----------
    path : str ::
        Path of the file provided as string.

    RETURNS
    -------
    numpy.array ::
        Numpy array file for the constants file.
    """
    return np.loadtxt(path, comments="%")


def load_client_data(path):
    """Utility function to load client-data files

    Utility funciton for loading files utilising the
    numpy binary file as storage format for the data
    for simulation.

    PARAMETERS
    ----------
    path : str ::
        Path of the file provided as string.

    RETURNS
    -------
    numpy.array ::
        Numpy array file for the constants file.
    """
    return np.load(path)


def load_data(path):
    """Utility Function to load data components

    Utility function for loading the required constants
    and feed data from the specified path. The function
    automating chooses between the plaintext parser and
    the numpy binary format parser.

    PARAMETERS
    ----------
    path : str ::
        Path of the file provided as string.

    RETURNS
    -------
    numpy.array ::
        Numpy array file for the constants file.
    """
    try:
        return np.loadtxt(path, comments="%")
    except UnicodeDecodeError as e:
        return np.load(path)
