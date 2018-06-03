#!/usr/bin/env/python

import numpy as np

"""
Utility

A collection of utility functions that can
be used by multiple scripts if needed. This is just for
allowing different tools accesses to common functions.
"""


def load_constants(path):
    """Utility function to load constants files"""
    return np.loadtxt(path, comments="%")


def load_client_data(path):
    """Utility function to load client-data files"""
    return np.load(path)

