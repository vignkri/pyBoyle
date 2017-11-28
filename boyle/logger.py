#!/usr/bin/env python

import os
import logging

"""
Data Logging

Toolset for logging Simulation and Experiment Setup
"""


home = os.getenv("HOME")
logging_folder = ".BoyleLogs"
try:
    if not os.path.exists(os.path.join(home, logging_folder)):
        os.mkdir(os.path.join(home, logging_folder))
    assert os.path.exists(os.path.join(home, logging_folder))
    path = os.path.join(home, logging_folder)
except OSError as e:
    raise

fh = logging.FileHandler(os.path.join(path, "Boyle.log"))
_format = logging.Formatter("%(asctime)s-%(levelname)s-%(module)s-%(message)s")
fh.setLevel(logging.DEBUG)
fh.setFormatter(_format)

simulationLogger = logging.getLogger(__name__)
simulationLogger.setLevel(logging.DEBUG)
simulationLogger.addHandler(fh)
