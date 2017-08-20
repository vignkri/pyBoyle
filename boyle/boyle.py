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

# Define reaction rates and growth factors
k0_carbon = dataset.mu_max.get("params").get("k0_carbon")
k0_prot = dataset.mu_max.get("params").get("k0_prot")
# --
mu_max = dataset.mu_max.get("params").get("mu_max")
mu_max_t0 = dataset.mu_max.get("params").get("mu_max_t0")

# --
boyle_logger.info("Finished setting up constants.")

# Constant One Argument
names = ["ks", "ks_nh3", "pk_low", "pk_high", "ks_nh3", "ki_carbon",
         "ki_prot", "ki_hac_hpr", "ki_hac_hbut", "ki_hac_hval",
         "ki_nh3_hac", "ki_lcfa"]
constants_one = [dataset.Const1.get("params").get(item) for item in names]

# XX-Val Argument
names = ["k_h", "ka_nh4", "ka_hac", "ka_hpr", "ka_hbut", "ka_hval",
         "ka1_co2", "ka2_co2", "ka_h2s", "ka_h2po4", "kw"]
xxval = [dataset.henry_constants.get(item) for item in names]

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
