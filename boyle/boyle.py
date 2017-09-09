#!/usr/bin/python

import yaml

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


# Regulate Settings
temp = regulate_settings.get("temperature")
flow_in = regulate_settings.get("flow_in")
flow_out = regulate_settings.get("flow_out")

# Days to hourly flow rate change
flow_in = flow_in / 24
flow_out = flow_out / 24

# Import datasets
dataset = Parameters(configuration)
dataset.process_data(temp=temp, flow_in=flow_in, flow_out=flow_out)
boyle_logger.info("Input data loaded.")

# --
boyle_logger.info("Finished setting up constants.")

_config = dict(initial=dataset.Initial.get("value"),
               start_time=start_time, end_time=end_time,
               step=step_size, metadata=experiment_name)

solver = Manager(model.standard, config=_config)
solver.initialize_solver(iname="vode", i_params=dataset._solver)
solver.function_parameters(parameters=[dataset])
solver.start()
