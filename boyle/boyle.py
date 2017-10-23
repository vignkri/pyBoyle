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


# Import datasets
dataset = Parameters(configuration)
boyle_logger.info("Input data loaded.")

# --
boyle_logger.info("Finished setting up constants.")

solver = Manager(model.standard, frame=dataset)
solver.start()
