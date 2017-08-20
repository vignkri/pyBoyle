#!/usr/bin/python

import yaml
import numpy as np

import model
from frame import Parameters
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
solver_settings = configuration.get("solver")
regulate_settings = configuration.get("regulate")
boyle_logger.info("Configuration file loaded.")

# Simulation settings
experiment_name = configuration.get("name")
start_time = settings.get("t_initial")
end_time = settings.get("t_final")
step_size = settings.get("step_size")

# Solver Settings
solver_method = solver_settings.get("method")
solver_order = solver_settings.get("order")
solver_nsteps = solver_settings.get("nsteps")
absolute_tolerance = solver_settings.get("absolute")
relative_tolerance = solver_settings.get("relative")

# Regulate Settings
temp = regulate_settings.get("temperature")
flow_in = regulate_settings.get("flow_in")
flow_out = regulate_settings.get("flow_out")

# Import datasets
dataset = Parameters("./sample")
const1 = np.loadtxt("./sample/Const1", comments="%")
const2 = np.loadtxt("./sample/Const2", comments="%")
boyle_logger.info("Input data loaded.")

substrate_inflow = flow_in * dataset.regulate.get("value")[4:]
dataset.process_data(temp=temp, flow_in=flow_in, flow_out=flow_out)

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
k0_carbon = dataset.mu_max.get("params").get("k0_carbon")
k0_prot = dataset.mu_max.get("params").get("k0_prot")
# --
mu_max = dataset.mu_max.get("params").get("mu_max")
mu_max_t0 = dataset.mu_max.get("params").get("mu_max_t0")

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
# XX-Val Argument
xxval = [k_h, ka_nh4, ka_hac, ka_hpr, ka_hbut, ka_hval, ka1_co2,
         ka2_co2, ka_h2s, ka_h2po4, kw]

_config = dict(initial=dataset.Initial.get("value"),
               start_time=start_time, end_time=end_time,
               step=step_size, metadata=experiment_name)
_solver_params = dict(method=solver_method, order=solver_order,
                      rtol=relative_tolerance, atol=absolute_tolerance,
                      nsteps=solver_nsteps)
_parameters = [constants_one, mu_max, xxval, mu_max_t0,
               [k0_carbon, k0_prot], [flow_in, flow_out],
               dataset.yc.get("value"), substrate_inflow]

solver = Manager(model.standard, config=_config)
solver.initialize_solver(iname="vode", i_params=_solver_params)
solver.function_parameters(parameters=_parameters)
solver.start()
