#!/usr/bin/env python

import logging

"""
Data Logging

Toolset for logging Simulation and Experiment Setup
"""

boyle_logger = logging.getLogger("boyle")
boyle_logger.setLevel(logging.DEBUG)

manager_logger = logging.getLogger("manager")
manager_logger.setLevel(logging.DEBUG)

model_logger = logging.getLogger("manager")
model_logger.setLevel(logging.DEBUG)

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
# --
boyle_logger.addHandler(fh)
boyle_logger.addHandler(ch)
# --
manager_logger.addHandler(fh)
manager_logger.addHandler(ch)
# --
model_logger.addHandler(fh)
model_logger.addHandler(ch)
