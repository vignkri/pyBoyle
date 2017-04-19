#!/usr/bin/python

import csv
import numpy as np

"""
Standard Computation Model
"""


def standard(time, y0,
             logging_headers, constant_ones, mu_max, xxval, mu_max_t0,
             k0_zeros, flow, yieldc, inflow):
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
    volume = y0[0]
    degraders = y0[21:29]
    substrate = y0[1:20]

    # -- LOGGER --
    subheader = logging_headers.get("substrate")
    subvalues = dict(zip(subheader, [time] + y0[0:20].tolist()))
    with open("./logging/substrates.log", "a") as sublog:
        sublogwriter = csv.DictWriter(sublog, fieldnames=subheader)
        sublogwriter.writerow(subvalues)
    # -- LOGGER --
    degrheader = logging_headers.get("degrader")
    degrvalues = dict(zip(degrheader, [time] + degraders.tolist()))
    with open("./logging/degraders.log", "a") as degrlog:
        degrlogwriter = csv.DictWriter(degrlog, fieldnames=degrheader)
        degrlogwriter.writerow(degrvalues)
    # -- LOGGER ALL VALUES --
    allheader = subheader + degrheader[1:] + ["g1", "g2", "g3", "g4"]
    subvalues = dict(zip(allheader, [time] + y0.tolist()))
    with open("./logging/all.log", "a") as initlog:
        initwriter = csv.DictWriter(initlog, fieldnames=allheader)
        initwriter.writerow(subvalues)

    # y0 Chemical Concentrations
    carbo_is, carbo_in, carbon, lipids, lcfa, \
        prot_is, prot_in, amino, nh3, hac, hpr, hbut, hval, \
        ch4, co2, h2s, z, h2po4, a = substrate
    dead_cells = y0[20]

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
    Hfunc = 1
    # --
    while abs(Hfunc - H) > 1e-12:
        Hfunc = co2 / 44 * ka1_co2 * (H + 2 * ka2_co2) / \
            (H * (H + ka1_co2) + ka1_co2 * ka2_co2) + \
            hac / 60 * ka_hac / (H + ka_hac) + \
            hpr / 74 * ka_hpr / (H + ka_hpr) + \
            hbut / 88 * ka_hbut / (H + ka_hbut) + \
            hval / 102 * ka_hval / (H + ka_hval) + \
            a / 35.5 + \
            h2po4 / 31 * (H + 2 * ka_h2po4) / (H + ka_h2po4) - \
            nh3 / 14 * H / (H + ka_nh4) - \
            z / 39 + \
            kw / H
        # --
        dhfunc_dh = - co2 / 44 * ka1_co2 * (H * (H + 4 * ka2_co2) + ka1_co2 *
                                            ka2_co2) / \
            (H * (H + ka1_co2) + ka1_co2 * ka2_co2)**2 - \
            hac / 60 * ka_hac / (H + ka_hac)**2 - \
            hpr / 74 * ka_hpr / (H + ka_hpr)**2 - \
            hbut / 88 * ka_hbut / (H + ka_hbut)**2 - \
            hval / 102 * ka_hval / (H + ka_hval)**2 - \
            kw / (H * H) - \
            h2po4 / 31 * ka_h2po4 / (H + ka_h2po4)**2 - \
            nh3 / 14 * ka_nh4 / (H + ka_nh4)**2
        H = H - (Hfunc - H) / (dhfunc_dh - 1)
    # --
    pH = -np.log10(H)

    # Calculation of growth rates
    f_ph = (1 + 2 * 10**(0.5 * (pk_low - pk_high))) / \
        (1 + 10**(pH - pk_high) + 10**(pk_low - pH))
    f_ph = f_ph.reshape(-1, 1)

    mu = mu_max * f_ph

    # -- LOGGING --
    muheader = logging_headers.get("mu")
    muvalues = dict(zip(muheader, [time] + mu[:, 0].tolist() + [False]))
    with open("./logging/mu_values.log", "a") as muvallog:
        muwriter = csv.DictWriter(muvallog, fieldnames=muheader)
        muwriter.writerow(muvalues)

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
    mu[6, 0] = mu[6, 0] * hval * nh3 * ki_lcfa[6] * ki_hac_hval / \
        ((ks[6] + hval) * (ks_nh3[6] + nh3) * (lcfa + ki_lcfa[6]) *
         (hac + ki_hac_hval))
    #
    mu[7, 0] = mu[7, 0] * hac * nh3 * ki_lcfa[7] * ki_nh3_hac / \
        ((ks[7] + hac) * (ks_nh3[7] + nh3) * (lcfa + ki_lcfa[7]) *
         (nh3 * ka_nh4 / (H + ka_nh4) + ki_nh3_hac))

    # -- LOGGING --
    muheader = logging_headers.get("mu")
    muvalues = dict(zip(muheader, [time] + mu[:, 0].tolist() + [True]))
    with open("./logging/mu_values.log", "a") as muvallog:
        muwriter = csv.DictWriter(muvallog, fieldnames=muheader)
        muwriter.writerow(muvalues)

    # Calculate growth, death and reaction rates
    cell_death = (mu_max_t0 * degraders.reshape(-1, 1)) * kd0
    cell_decay = 0.01 * dead_cells
    # -- column of z
    z = np.array([
        cell_decay,
        carbo_is * k0_carbon * ki_carbon / (ki_carbon + hac + 0.811 * hpr +
                                            0.659 * hbut),
        prot_is * k0_prot * ki_prot / (ki_prot + hac + 0.811 * hpr +
                                       0.659 * hbut)
    ])
    z_two = (mu * degraders.reshape(-1, 1)).reshape(-1)
    z = np.concatenate((z, z_two))

    # Computation of Gasflow
    flow_in, flow_out = flow
    y_dot = np.zeros((33,))
    y_dot[0] = flow_in - flow_out
    y_dot[1:17] = (yieldc.conj().transpose()).dot(z)
    y_dot[20] = np.sum(cell_death) - cell_decay
    y_dot[21:29] = z[3:] - cell_death.reshape(-1)
    y_dot[1:29] = y_dot[1:29] + (inflow - flow_in * y0[1:29]) / volume

    # Calculation of gasflow
    molar_mass = np.array([14, 16, 44, 34])
    conc = np.array([nh3, ch4, co2, h2s]) / molar_mass
    dconc_dt = np.array([y_dot[[9, 14, 15, 16]]]) / molar_mass
    # --
    a = np.array([ka_nh4 / (H + ka_nh4),
                  1,
                  H * H / (H * (H + ka1_co2) + ka1_co2 * ka2_co2),
                  H / (H + ka_h2s)]) / k_h
    da_dH = np.array([
        -ka_nh4 / (H + ka_nh4)**2,
        0,
        ka1_co2 * H * (H + 2 * ka2_co2) / (H * (H + ka1_co2) +
                                           ka1_co2 * ka2_co2)**2,
        ka_h2s / (H + ka_h2s)**2
    ]) / k_h
    # --
    dH_dt = - (ka_hac / (ka_hac + H) * y_dot[9] / 60 +
               ka_hpr / (ka_hac + H) * y_dot[10] / 74 +
               ka_hbut / (ka_hac + H) * y_dot[11] / 88 +
               ka_hval / (ka_hac + H) * y_dot[12] / 102 +
               y_dot[18] / 35.5 -
               y_dot[16] / 39 +
               (1 + ka_h2po4 / (ka_h2po4 - H)) * y_dot[17] / 31) / \
        (
            ((ka1_co2 - 1) * ka2_co2 - H * H) * co2 / 44 /
            (H * (H + ka1_co2) + ka1_co2 * ka2_co2)**2 +
            ka_hac / (ka_hac + H)**2 * hac / 60 +
            ka_hpr / (ka_hac + H)**2 * hpr / 74 +
            ka_hbut / (ka_hac + H)**2 * hbut / 88 +
            ka_hval / (ka_hac + H)**2 * hval / 102 -
            kw / (H * H) - 1 +
            ka_h2po4 / (ka_h2po4 + H)**2 * h2po4 / 31 +
            ka_nh4 / (ka_nh4 + H)**2 * nh3 / 14)
    # --
    gasflow = np.sum(a * dconc_dt + da_dH * dH_dt * conc) / \
        np.dot(np.sum(a * a * conc), a) * conc
    gasloss = gasflow * molar_mass
    gasflow = (gasflow.dot(volume)).dot(22.4)
    # --
    y_dot[-4:] = gasflow
    y_dot[9] = y_dot[9] - gasloss[0]
    y_dot[[14, 15, 16]] = y_dot[[14, 15, 16]] - gasloss[[1, 2, 3]]
    # --
    y_dot = y_dot.reshape(-1)

    # -- LOGGING --
    dydtheader = logging_headers.get("dydt")
    dydtvalues = dict(zip(dydtheader, [time] + y_dot.tolist()))
    with open("./logging/dydt.log", "a") as dylog:
        dylogwriter = csv.DictWriter(dylog, fieldnames=dydtheader)
        dylogwriter.writerow(dydtvalues)
    # --

    return y_dot
