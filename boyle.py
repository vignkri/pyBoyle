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

# Set up regulation
start_time = 0
end_time = regulate[0]
step_size = 0.5

# Substrate Conditions
temp = regulate[1]
flow_in = regulate[2]
flow_out = regulate[3]
substrate_inflow = flow_in * regulate[4:]

# Compute Temperature Dependent Constants
mu_max = np.zeros((10, 1))
mu_max_t0 = np.zeros((10, 1))
for index in range(0, 10):
    mu_max_t0[index] = const1[index, 0]
    alpha = const1[index, 1]
    t0 = const1[index, 2]
    t_opt = const1[index, 3]
    t_max = const1[index, 4]
    #
    if temp < t_opt:
        mu_max[index] = mu_max_t0[index] + alpha * (temp - t0)
    else:
        mu_max[index] = (mu_max_t0[index] + alpha * (t_opt - t0)) * \
            (t_max - temp) / (t_max - t_opt)
# --
