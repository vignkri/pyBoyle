#!/usr/bin/python

import numpy as np


"""
Boyle Python

Simulation agent for Biogas Production
"""


# Import datasets
initial = np.loadtxt("./sample/Initial", comments="%")
yield_c = np.loadtxt("./sample/yc", comments="%")
regulate = np.loadtxt("./sample/regulate", comments="%")
const1 = np.loadtxt("./sample/Const1", comments="%")
const2 = np.loadtxt("./sample/Const2", comments="%")

# Set up initial values
volume = initial[0]
substrate = initial[1:20]
degraders = initial[20:29]
gas_conc = np.zeros((4,))
