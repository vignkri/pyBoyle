#!/usr/bin/python

import csv
import numpy as np

"""
Standard Computation Model
"""


def standard(initial, time,
             logging_headers, constant_ones, mu_max):
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
    # -- LOGGER --
    subheader = logging_headers.get("substrate")
    subvalues = dict(zip(subheader, initial[0:20]))
    with open("./logging/substrates.log", "a") as sublog:
        sublogwriter = csv.DictWriter(sublog, fieldnames=subheader)
        sublogwriter.writerow(subvalues)
    # -- LOGGER --
    degrheader = logging_headers.get("degrader")
    degrvalues = dict(zip(degrheader, degraders))
    with open("./logging/degraders.log", "a") as degrlog:
        degrlogwriter = csv.DictWriter(degrlog, fieldnames=degrheader)
        degrlogwriter.writerow(degrvalues)

    # Initial Chemical Concentrations
    carbo_is, carbo_in, carbon, lipids, lcfa, \
        prot_is, prot_in, amino, nh3, hac, hpr, hbut, hval, \
        ch4, co2, h2s, z, h2po4, a = substrate
    dead_cells = degraders[0]
    degraders = degraders[1:]

    # Constant One Argument
    ks, ks_nh3, pk_low, pk_high, ks_nh3, \
        ki_carbon, ki_prot, ki_hac_hpr, ki_hac_hbut, \
        ki_hac_hval, ki_nh3_hac, ki_lcfa = constant_ones

    # Set up pH Computation
    H = 1e-8
    pH = -np.log10(H)

    # Calculation of growth rates
    f_ph = (1 + 2 * 10**(0.5 * (pk_low - pk_high))) / \
        (1 + 10**(pH - pk_high) + 10**(pk_low - pH))
    mu = mu_max * np.atleast_2d(f_ph).T
