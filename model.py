#!/usr/bin/python

import numpy as np

"""
Standard Computation Model
"""


def standard(initial, time,
             logging_headers):
    """Standard Integrator Model

    PARAMETERS
    ----------
    init : list
        A list of four values:
            - volume
            - substrates : numpy.array
            - degraders  : numpy.array
            - gas_values : numpy.array

    OTHER PARAMETERS
    ----------------
    """
    # -- set up parts of values
    volume = initial[0]
    substrate = initial[1:20]
    degraders = initial[20:29]
    gas_conc = initial[-4:]
