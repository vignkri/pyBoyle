#!/usr/bin/python

import csv

"""
Simulation Logging Framework

A custom framework to handle the creation
of log files, data files during simulation for
further analysis.
"""


class Simulog:
    def __init__(self, path):
        """Initialize data logging folder path."""
        self._path = path
