#!/usr/bin/env/python


"""
pH computation

A series of functions dedicated to
computing the pH value of the input feed
from the model. The model calls the functions
in this module in order to compute the pH
for the feed and the inoculum present in
the model at that particular time.
"""


def newton_raphson(H, Hfunc, i=0, **kwargs):
    """pH Computation using Newton-Raphson method"""
    # -- set up variables
    co2, ka1_co2, ka2_co2 = kwargs.get("co2")
    hac, ka_hac = kwargs.get("HAc")
    hpr, ka_hpr = kwargs.get("HPr")
    hbut, ka_hbut = kwargs.get("HBut")
    hval, ka_hval = kwargs.get("HVal")
    a, z, kw = kwargs.get("Other")
    h2po4, ka_h2po4 = kwargs.get("h2po4")
    nh3, ka_nh4 = kwargs.get("NH3")
    # --
    Hfunc = co2 / 44 * ka1_co2 * (H + 2 * ka2_co2) / (H * (H + ka1_co2) +
                                                      ka1_co2 * ka2_co2) + \
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
                                        ka2_co2) / (H * (H + ka1_co2) +
                                                    ka1_co2 * ka2_co2)**2 - \
        hac / 60 * ka_hac / (H + ka_hac)**2 - \
        hpr / 74 * ka_hpr / (H + ka_hpr)**2 - \
        hbut / 88 * ka_hbut / (H + ka_hbut)**2 - \
        hval / 102 * ka_hval / (H + ka_hval)**2 - \
        kw / (H**2) - \
        h2po4 / 31 * ka_h2po4 / (H + ka_h2po4)**2 - \
        nh3 / 14 * ka_nh4 / (H + ka_nh4)**2
    # --
    H = H - (Hfunc - H) / (dhfunc_dh - 1)
    return H, Hfunc
