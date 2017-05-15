#!/usr/bin/env python

import logging

"""
Data Logging

Toolset for logging Simulation and Experiment Setup
"""

logger = logging.getLogger("boyle")
logger.setLevel(logging.DEBUG)
# --
fh = logging.FileHandler("./logs/simulation.log")
fh.setLevel(logging.INFO)
# --
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter
formatter = logging.Formatter('%(asctime)s-%(levelname)s-%(name)s %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)
