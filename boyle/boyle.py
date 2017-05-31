#!/usr/bin/python

import numpy as np
import scipy.integrate

import model
import export
from scribe import logger

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
logger.info("Input data loaded.")

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

# Set up Const1 Parameters
kd0 = 0.05
ks = const1[2:, 5]
ks_nh3 = const1[2:, 6]
pk_low = const1[2:, 9]
pk_high = const1[2:, 10]
# --
ki_carbon = const1[0, 7]
ki_prot = const1[1, 7]
ki_hac_hpr = const1[6, 7]
ki_hac_hbut = const1[7, 7]
ki_hac_hval = const1[8, 7]
ki_nh3_hac = const1[9, 8]
ki_lcfa = const1[2:, 8]

# Define reaction rates and growth factors
k0_carbon = mu_max[0, 0]
k0_prot = mu_max[1, 0]
# --
mu_max_t0 = mu_max_t0[2:]
mu_max = mu_max[2:]

# Defining Henry Constant Values
delta_temp = temp - const2[:, 1]
henry_constants = const2[:, 0] + delta_temp * const2[:, 2] + \
    delta_temp**2 * const2[:, 3] + delta_temp**3 * const2[:, 4]
# --
k_h = henry_constants[[1, 7, 8, 11]]
# -- log inverse values
ka1_lcfa = 10**(-henry_constants[0])
ka_nh4 = 10**(-henry_constants[2])
ka_hac = 10**(-henry_constants[3])
ka_hpr = 10**(-henry_constants[4])
ka_hbut = 10**(-henry_constants[5])
ka_hval = 10**(-henry_constants[6])
ka1_co2 = 10**(-henry_constants[9])
ka2_co2 = 10**(-henry_constants[10])
ka_h2s = 10**(-henry_constants[12])
ka_h2po4 = 10**(-henry_constants[13])
kw = 10**(-henry_constants[14])
# --
logger.info("Finished setting up constants.")

# Create Logging Parameters
simport = export.Simulation()
logger.info("Create simulation data exporter.")

# Constant One Argument
constants_one = [ks, ks_nh3, pk_low, pk_high, ks_nh3,
                 ki_carbon, ki_prot, ki_hac_hpr, ki_hac_hbut,
                 ki_hac_hval, ki_nh3_hac, ki_lcfa]
# XXVal Argument
xxval = [k_h, ka_nh4, ka_hac, ka_hpr, ka_hbut, ka_hval, ka1_co2,
         ka2_co2, ka_h2s, ka_h2po4, kw]
# Set up integrator
time_array = np.linspace(start_time, end_time,
                         (end_time - start_time)/step_size)
initial_values = np.concatenate((np.array([volume]), substrate,
                                 degraders, gas_conc))

solver = scipy.integrate.ode(model.standard) \
    .set_integrator("vode", method="bdf", order=1, rtol=1e-4, atol=1e-8,
                    nsteps=2)
solver.set_initial_value(y=initial_values, t=start_time)
solver.set_f_params(constants_one, mu_max, xxval, mu_max_t0,
                    [k0_carbon, k0_prot], [flow_in, flow_out], yield_c,
                    substrate_inflow, simport)
logger.info("Set up integrator.")
logger.info("Starting Integration.")

result_set = []
while solver.successful() and solver.t < end_time:
    logger.info("Computation Time %.2f" % solver.t)
    y_dot = solver.integrate(solver.t + step_size, step=True)
    simport._append_values("result", [solver.t] + list(y_dot))
    result_set.append(np.append(np.array([solver.t]), (y_dot)))

for idx in reversed(list(range(1, len(result_set)))):
    result_set[idx][29:] = (result_set[idx][29:] - result_set[idx-1][29:]) / (
        result_set[idx][0] - result_set[idx-1][0]
    ) / 1000
    simport._append_values("processed", result_set[idx])
