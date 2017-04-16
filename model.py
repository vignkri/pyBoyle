#!/usr/bin/python

import csv
import numpy as np

"""
Standard Computation Model
"""


def standard(initial, time,
             logging_headers, constant_ones):
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
    with open("./logging/substrate.log", "a") as sublog:
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
