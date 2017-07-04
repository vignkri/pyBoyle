#!/usr/bin/python

import os
import csv
import time

"""
Simulation Logging Framework

A custom framework to handle the creation
of log files, data files during simulation for
further analysis.
"""


class BoyleOutput(object):
    """Resultant export object from the simulation."""
    def __init__(self, exp_name):
        self._name = exp_name
        self._run_date = time.gmtime()
        # result key_list
        self.__available_headers = dict(
            mu=["time", "one", "two", "three", "four", "five", "six", "seven",
                "eight", "PostFlag"],
            substrates=["time", "volume", "carbo_is", "carbo_in", "carbon",
                        "lipids", "lcfa", "prot_is", "prot_in", "amino",
                        "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                        "co2", "h2s", "zplus", "h2po4", "aminus"],
            degraders=["time", "carb_degr", "amino_degr", "lipid_degr",
                       "lcfa_degr", "prop_degr", "butyr_degr",
                       "valer_degr", "acet_degr"],
            result=["time", "volume", "carbois", "carboin", "carbon",
                    "lipids", "lcfa", "protis", "protin", "amino",
                    "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                    "co2", "h2s", "zplus", "h2po4", "aminus", "deadcell",
                    "carb_degr", "amino_degr", "lipid_degr", "lcfa_degr",
                    "prop_degr", "butyr_degr", "valer_degr", "acet_degr",
                    "gfnh3", "gfch4", "gfco2", "gfh2s"],
            ph=["time", "ph"]
        )

    def _update(self, attrib, value):
        """Update the attribute if the value exists"""
        try:
            getattr(self, attrib)
        except AttributeError as e:
            setattr(self, attrib, [self.__available_headers.get(attrib)])
            # -- log the error in the logging file.
        else:
            getattr(self, attrib).append(value)


class Export:
    def __init__(self, meta):
        """Initialize data logging folder path."""
        try:
            base_path = "./logs"
            if not os.path.exists(base_path):
                os.mkdir(base_path)
            self._path = os.path.join(base_path, meta)
            if not os.path.exists(self._path):
                os.mkdir(os.path.join(base_path, meta))
        except:
            raise
        # Log date as folder
        self.itime = time.gmtime()
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
                       "gfnh3", "gfch4", "gfco2", "gfh2s"],
            "ph": ["time", "ph"],
            "processed": ["time", "volume", "carbois", "carboin", "carbon",
                          "lipids", "lcfa", "protis", "protin", "amino",
                          "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                          "co2", "h2s", "zplus", "h2po4", "aminus", "deadcell",
                          "carb_degr", "amino_degr", "lipid_degr", "lcfa_degr",
                          "prop_degr", "butyr_degr", "valer_degr", "acet_degr",
                          "gfnh3", "gfch4", "gfco2", "gfh2s"],
        }
        # Create file names
        self.__file_names = {
            "mu": os.path.join(self._path, "mu.dat"),
            "substrates": os.path.join(self._path, "substrates.dat"),
            "degraders": os.path.join(self._path, "degraders.dat"),
            "result": os.path.join(self._path, "result.dat"),
            "ph": os.path.join(self._path, "ph.dat"),
            "processed": os.path.join(self._path, "processed.dat")
        }
        # Initialize creation of files
        self.__create_files()

    def __create_files(self):
        """Create files with custom headers."""
        for key in self.__file_names.keys():
            file_path = self.__file_names.get(key)
            fnames = self.__headers.get(key)
            try:
                with open(file_path, "w") as open_log_file:
                    dwriter = csv.DictWriter(open_log_file, fieldnames=fnames)
                    dwriter.writeheader()
            except:
                raise

    def _append_values(self, sim_locn, values):
        """Append values to the files from the dataset"""
        _header = self.__headers.get(sim_locn)
        file_path = self.__file_names.get(sim_locn)
        a_values = dict(zip(_header, values))
        try:
            with open(file_path, "a") as open_log_file:
                dwriter = csv.DictWriter(open_log_file, fieldnames=_header)
                dwriter.writerow(a_values)
        except:
            raise
