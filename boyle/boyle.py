#!/usr/bin/python

import os
import yaml
import argparse

import model
from frame import Parameters
from logger import simulationLogger
from manager import Manager

"""
Boyle Python

Simulation agent for Biogas Production
"""


def main(cfg_path):
    config_file = cfg_path
    with open(config_file, "r") as config_stream:
        configuration = yaml.load(config_stream)
    # --
    # Import datasets
    dataset = Parameters(configuration)
    simulationLogger.info("Loaded Input data into frame.")

    solver = Manager(model.standard, frame=dataset)
    solver.start()


if __name__ == "__main__":
    cwd = os.getcwd()
    # -- Set up argument parser
    parser = argparse.ArgumentParser(description="Simulation Engine")
    parser.add_argument("-path", type=str,
                        help="Path of simulation configuration file.")
    try:
        args = parser.parse_args()
        assert args.path is not None
    except AssertionError as e:
        print("Please provide path of configuration file.")
        simulationLogger.critical("No arguments are provided. %s" % e)
        raise
    # --
    try:
        assert os.path.isfile(os.path.join(cwd, args.path))
    except AssertionError as e:
        print("Configuration file not found.")
        simulationLogger.critical("Configuration file not found. %s" % e)
        raise
    finally:
        main(os.path.join(cwd, args.path))
