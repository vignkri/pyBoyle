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
        # Create internally useful headers
        self.__headers = {
            "mu": ["time", "mu_one", "mu_two", "mu_three", "mu_four",
                   "mu_five", "mu_six", "mu_seven", "mu_eight", "afterflag"],
            "substrates": ["time", "volume", "carbo_is", "carbo_in", "carbon",
                           "lipids", "lcfa", "prot_is", "prot_in", "amino",
                           "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                           "co2", "h2s", "zplus", "h2po4", "aminus"],
            "degraders": ["time", "carb_degr", "amino_degr", "lipid_degr",
                          "lcfa_degr", "prop_degr", "butyr_degr",
                          "valer_degr", "acet_degr"],
            "result": ["time", "volume", "carbois", "carboin", "carbon",
                       "lipids", "lcfa", "protis", "protin", "amino",
                       "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                       "co2", "h2s", "zplus", "h2po4", "aminus", "deadcell",
                       "carb_degr", "amino_degr", "lipid_degr", "lcfa_degr",
                       "prop_degr", "butyr_degr", "valer_degr", "acet_degr",
                       "gfnh3", "gfch4", "gfco2", "gfh2s"]
        }

    def _create_files(self, path, fnames):
        """Create files with custom headers."""
        try:
            with open(path, "w") as open_log_file:
                dwriter = csv.DictWriter(open_log_file, fieldnames=fnames)
                dwriter.writeheader()
        except:
            raise
