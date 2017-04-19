#!/usr/bin/python

import csv
import numpy as np
import pandas as pd
import scipy.integrate

import model

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

# Create Logging Parameters
# Logging for mu
mu_header = ["time", "mu_one", "mu_two", "mu_three", "mu_four",
             "mu_five", "mu_six", "mu_seven", "mu_eight", "after_flag"]
with open("./logging/mu_values.log", "w") as mu_val_log:
    mu_writer = csv.DictWriter(mu_val_log, fieldnames=mu_header)
    mu_writer.writeheader()
# Logging for Substrate Values
substrate_header = ["time", "volume",
                    "carbo_is", "carbo_in", "carbon",
                    "lipids", "lcfa", "prot_is", "prot_in", "amino",
                    "nh3", "hac", "hpr", "hbut", "hval", "ch4", "co2",
                    "h2s", "zplus", "h2po4", "aminus"]
with open("./logging/substrates.log", "w") as substrate_log:
    subs_writer = csv.DictWriter(substrate_log, fieldnames=substrate_header)
    subs_writer.writeheader()
# Logging for Degrader Values
degrader_header = ["time", "carb_degr", "amino_degr", "lipid_degr",
                   "lcfa_degr", "prop_degr", "butyr_degr", "valer_degr",
                   "acet_degr"]
with open("./logging/degraders.log", "w") as degrader_log:
    degs_writer = csv.DictWriter(degrader_log, fieldnames=degrader_header)
    degs_writer.writeheader()
# Logging for All Values
all_header = substrate_header + ["deadcell"] + degrader_header[1:] + \
    ["g1", "g2", "g3", "g4"]
with open("./logging/all.log", "w") as all_log:
    all_writer = csv.DictWriter(all_log, fieldnames=all_header)
    all_writer.writeheader()
# Loggin for Change values
dydt_header = substrate_header + degrader_header[1:]
with open("./logging/dydt.log", "w") as dydt_log:
    dydt_writer = csv.DictWriter(dydt_log, fieldnames=dydt_header)
    dydt_writer.writeheader()

# Logging Headers Argument
loggers = {"substrate": substrate_header,
           "degrader": degrader_header,
           "dydt": dydt_header,
           "mu": mu_header}
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
solver.set_f_params(loggers, constants_one, mu_max, xxval, mu_max_t0,
                    [k0_carbon, k0_prot], [flow_in, flow_out], yield_c,
                    substrate_inflow)

# -- List of outputs
results_ = []
while solver.successful() and solver.t < end_time:
    print("Computation Time: %.2f" % solver.t)
    y_dot = solver.integrate(solver.t + step_size, step=True)
    result = dict(zip(all_header, np.append(solver.t, y_dot)))
    results_.append(result)

resultFrame = pd.DataFrame.from_records(results_)
resultFrame.to_csv("./logging/results.log")
