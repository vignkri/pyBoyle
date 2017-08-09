#!/usr/bin/python

import yaml
import numpy as np

import model
from logger import boyle_logger
from manager import Manager

"""
Boyle Python

Simulation agent for Biogas Production
"""

# Configuration Load
config_file = "./simulation.yaml"
with open(config_file, "r") as config_stream:
    configuration = yaml.load(config_stream)
# --
settings = configuration.get("settings")
# regulate = configuration.get("regulate")
boyle_logger.info("Configuration file loaded.")

# Simulation settings
start_time = settings.get("t_initial")
end_time = settings.get("t_final")
step_size = settings.get("step_size")

# Regulation Settings
# temp = regulate.get("temperature")
# flow_in = regulate.get("flow_in")
# flow_out = regulate.get("flow_out")
# substrate_inflow = flow_in * np.array(list(regulate.get("flow").values()))

# Import datasets
initial = np.loadtxt("./sample/Initial", comments="%")
yield_c = np.loadtxt("./sample/yc", comments="%")
const1 = np.loadtxt("./sample/Const1", comments="%")
const2 = np.loadtxt("./sample/Const2", comments="%")
regulate = np.loadtxt("./sample/regulate", comments="%")
boyle_logger.info("Input data loaded.")

temp = regulate[1]
flow_in = regulate[2]
flow_out = regulate[3]
substrate_inflow = flow_in * regulate[4:]

# Set up initial values
initial = np.concatenate((initial, np.zeros(4, )))

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
boyle_logger.info("Finished setting up constants.")

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

parameters = [constants_one, mu_max, xxval, mu_max_t0,
              [k0_carbon, k0_prot], [flow_in, flow_out], yield_c,
              substrate_inflow]

solver = Manager(model.standard,
                 config=dict(
                     initial=initial,
                     start_time=0,
                     end_time=1000,
                     step=step_size,
                     metadata="Test"
                 ))
solver.initialize_solver(iname="vode",
                         i_params=dict(
                             method="bdf",
                             order=1,
                             rtol=1e-4,
                             atol=1e-8,
                             nsteps=2
                         ))
solver.function_parameters(parameters=parameters)
solver.start()
