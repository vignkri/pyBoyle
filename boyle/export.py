#!/usr/bin/python

import os
import dill
import time

from logger import simulationLogger

"""
Simulation Logging Framework

A custom framework to handle the creation
of log files, data files during simulation for
further analysis.
"""


class BoyleOutput(object):
    """Resultant export object from the simulation."""
    def __init__(self, exp_name, model, output):
        self._name = exp_name
        self._model_name = model
        self._creation_date = time.gmtime()
        # result key_list
        self.__available_headers = dict(
            debug=["time", "mu_1", "mu_2", "mu_3", "mu_4", "mu_5",
                   "mu_6", "mu_7", "mu_8", "ph", "volume", "carbois",
                   "carboin", "carbon", "lipids", "lcfa", "protis",
                   "protin", "amino", "nh3", "hac", "hpr", "hbut",
                   "hval", "ch4", "co2", "h2s", "zplus", "h2po4",
                   "aminus", "deadcell", "carb_degr", "amino_degr",
                   "lipid_degr", "lcfa_degr", "prop_degr", "butyr_degr",
                   "valer_degr", "acet_degr", "gfnh3", "gfch4",
                   "gfco2", "gfh2s"],
            solution=["time", "volume", "carbois", "carboin", "carbon",
                      "lipids", "lcfa", "protis", "protin", "amino",
                      "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                      "co2", "h2s", "zplus", "h2po4", "aminus", "deadcell",
                      "carb_degr", "amino_degr", "lipid_degr", "lcfa_degr",
                      "prop_degr", "butyr_degr", "valer_degr", "acet_degr",
                      "gfnh3", "gfch4", "gfco2", "gfh2s", "gasrate"],
        )
        try:
            _base_path = output
            output_fldr = os.path.join(_base_path, "output")
            if not os.path.exists(output_fldr):
                os.mkdir(output_fldr)
            self._path = output_fldr
        except OSError as e:
            simulationLogger.error("OSError: BoyleOutput Creation Module.")
            raise

    def __repr__(self):
        """Create representation of the object for introspection."""
        return ("Simulation {name} performed on {DD}/{MON}/{YYYY}-{HH}:{MM}"
                .format(name=self._name, DD=self._creation_date.tm_mday,
                        MON=self._creation_date.tm_mon,
                        YYYY=self._creation_date.tm_year,
                        HH=self._creation_date.tm_hour,
                        MM=self._creation_date.tm_min))

    def as_pickle(self):
        """Store the data as a pickle."""
        self._finish_date = time.gmtime()
        try:
            with open(os.path.join(self._path, "output.pkl"), "wb") as _file:
                dill.dump(self, _file)
        except OSError as e:
            simulationLogger.error("OSError: Output pickle file not found."
                                   "Create output file/folder to continue.")
            raise

    def _update(self, attrname, value):
        """Update the attribute if the value exists"""
        try:
            getattr(self, attrname)
        except AttributeError as e:
            setattr(self, attrname, [self.__available_headers.get(attrname)])
            simulationLogger.warning("BoyleOutput: attribute '%s' not found"
                                     % (attrname))
            # -- log the error in the logging file.
        else:
            getattr(self, attrname).append(value)
