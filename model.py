#!/usr/bin/python

import csv
import numpy as np

"""
Standard Computation Model
"""


def standard(initial, time,
             logging_headers, constant_ones, mu_max, xxval, mu_max_t0,
             k0_zeros, flow, yieldc):
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
    # Constants
    kd0 = 0.05
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

    # XXVal Argument
    k_h, ka_nh4, ka_hac, ka_hpr, ka_hbut, ka_hval, ka1_co2, \
        ka2_co2, ka_h2s, ka_h2po4, kw = xxval

    # K_zeros
    k0_carbon, k0_prot = k0_zeros

    # Set up pH Computation
    H = 1e-8
    pH = -np.log10(H)

    # Calculation of growth rates
    f_ph = (1 + 2 * 10**(0.5 * (pk_low - pk_high))) / \
        (1 + 10**(pH - pk_high) + 10**(pk_low - pH))
    mu = mu_max * np.atleast_2d(f_ph).T
    # --
    mu[0, 0] = mu[0, 0] * carbon * nh3 * ki_lcfa[0] / \
        ((ks[0] + carbon) * (ks_nh3[0] + nh3) * (lcfa + ki_lcfa[0]))
    #
    mu[1, 0] = mu[1, 0] * amino * ki_lcfa[1] / \
        ((ks[1] + amino) * (lcfa + ki_lcfa[1]))
    #
    mu[2, 0] = mu[2, 0] * lipids * nh3 * ki_lcfa[2] / \
        ((ks[2] + lipids) * (ks_nh3[2] + nh3) * (lcfa + ki_lcfa[2]))
    #
    mu[3, 0] = mu[3, 0] * nh3 * lcfa / \
        ((lcfa + ks[3] + lcfa * lcfa / ki_lcfa[3]) * (ks_nh3[3] + nh3))
    #
    mu[4, 0] = mu[4, 0] * hpr * nh3 * ki_lcfa[4] * ki_hac_hpr / \
        ((ks[4] + hpr) * (ks_nh3[4] + nh3) * (lcfa + ki_lcfa[4]) *
         (hac + ki_hac_hpr))
    #
    mu[5, 0] = mu[5, 0] * hbut * nh3 * ki_lcfa[5] * ki_hac_hbut / \
        ((ks[5] + hbut) * (ks_nh3[5] + nh3) * (lcfa + ki_lcfa[5]) *
         (hac + ki_hac_hbut))
    #
    mu[6, 0] = mu[6, 0] * hval * nh3 * ki_lcfa[6] + ki_hac_hval / \
        ((ks[6] + hval) * (ks_nh3[6] + nh3) * (lcfa + ki_lcfa[6]) *
         (hac + ki_hac_hval))
    #
    mu[7, 0] = mu[7, 0] * hac * nh3 * ki_lcfa[7] * ki_nh3_hac / \
        ((ks[7] + hac) * (ks_nh3[7] + nh3) * (lcfa + ki_lcfa[7]) *
         (nh3 * ka_nh4 / (H + ka_nh4) + ki_nh3_hac))

    # -- LOGGING --
    muheader = logging_headers.get("mu")
    muvalues = dict(zip(muheader, mu.tolist()))
    with open("./logging/mu_values.log", "a") as muvallog:
        muwriter = csv.DictWriter(muvallog, fieldnames=muheader)
        muwriter.writerow(muvalues)

    # Calculate growth, death and reaction rates
    cell_death = mu_max_t0 * degraders * kd0
    cell_decay = 0.01 * dead_cells
    # -- column of z
    z = np.array([
        cell_decay,
        carbo_is / k0_carbon * ki_carbon / (ki_carbon + hac + 0.811 * hpr +
                                            0.659 * hbut),
        prot_is * k0_prot * ki_prot / (ki_prot + hac + 0.811 * hpr +
                                       0.659 * hbut)
    ])
    z_two = (mu * np.atleast_2d(degraders).T).reshape(-1)
    z = np.concatenate((z, z_two))

    # Computation of Gasflow
    flow_in, flow_out = flow
    y_dot = np.zeros((33, 1))
    y_dot[0] = flow_in - flow_out
    y_dot[1:17] = yieldc.transpose().dot(z)
