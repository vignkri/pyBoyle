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

regulate_settings = configuration.get("regulate")
boyle_logger.info("Configuration file loaded.")

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

solver = Manager(model.standard, config=dataset._simulation_config)
solver.initialize_solver(iname="vode", i_params=dataset._solver)
solver.function_parameters(parameters=[dataset])
solver.start()
