#!/usr/bin/python

import os
import csv
import time
from string import Template

"""
Simulation Logging Framework

A custom framework to handle the creation
of log files, data files during simulation for
further analysis.
"""


class Simulog:
    def __init__(self):
        """Initialize data logging folder path."""
        self._path = "./logs"
        if not os.path.exists(self._path):
            os.mkdir(self._path)
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
        # Create file names
        self.__file_names = {"mu": self._path + "/mu.dat",
                             "substrates": self._path + "/substrates.dat",
                             "degraders": self._path + "/degraders.dat",
                             "result": self._path + "/result.dat"}
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


class InternalLog:
    def __init__(self):
        """Create internal logs for operational statuses"""
        self._path = "./logs/"
        if not os.path.exists(self._path):
            os.mkdir(self._path)
        # --
        self.itime = time.localtime()
        self.experiment_directory = Template(
            "${path}/${year}${mon}${day}_${h}${m}.ilog").substitute(
                path=self._path, year=self.itime.tm_year,
                mon=self.itime.tm_mon, day=self.itime.tm_mday,
                h=self.itime.tm_hour, m=self.itime.tm_min)
        if not os.path.exists(self.experiment_directory):
            os.mkdir(self.experiment_directory)
