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
data_folder = configuration.get("data")
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

# Days to hourly flow rate change
flow_in = flow_in / 24
flow_out = flow_out / 24

# Import datasets
dataset = Parameters(data_folder)
dataset.process_data(temp=temp, flow_in=flow_in, flow_out=flow_out)
boyle_logger.info("Input data loaded.")

# --
boyle_logger.info("Finished setting up constants.")

_config = dict(initial=dataset.Initial.get("value"),
               start_time=start_time, end_time=end_time,
               step=step_size, metadata=experiment_name)
_solver_params = dict(method=solver_method, order=solver_order,
                      rtol=relative_tolerance, atol=absolute_tolerance,
                      nsteps=solver_nsteps)

solver = Manager(model.standard, config=_config)
solver.initialize_solver(iname="vode", i_params=_solver_params)
solver.function_parameters(parameters=[dataset])
solver.start()
