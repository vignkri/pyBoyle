#!/usr/bin/env python

"""
Formulae

A collection formulae
that do not have specific module
names for categorisation
"""


def computeHenryConstant(arr, temp):
    """Compute HenryConstants with provided array"""
    const2 = arr
    delta_temp = temp - const2["t0"]
    # Uses the following columns: X(T0), a, b, c
    # This does the equation hc = T + dt * a + dt^2 * b + dt^3 * c
    henry_constants = const2["xt0"] + delta_temp * const2["a"] + \
        delta_temp**2 * const2["b"] + delta_temp**3 * const2["c"]
    # --
    hc = dict(
        k_h=henry_constants[[5, 7, 8, 11]],  # K_H results
        # -- log inverse values
        ka1_lcfa=10**(-henry_constants[0]),  # LCFA
        ka_nh4=10**(-henry_constants[6]),  # NH4+
        ka_hac=10**(-henry_constants[1]),  # HAC
        ka_hpr=10**(-henry_constants[2]),  # HPr
        ka_hbut=10**(-henry_constants[3]),  # HBut
        ka_hval=10**(-henry_constants[4]),  # HVal
        ka1_co2=10**(-henry_constants[9]),  # pKa1 CO2
        ka2_co2=10**(-henry_constants[10]),  # pKa2 CO2
        ka_h2s=10**(-henry_constants[12]),  # pKa H2S
        ka_h2po4=10**(-henry_constants[13]),  # pKa H2PO4-
        kw=10**(-henry_constants[14])
    )
    return (hc, henry_constants)
