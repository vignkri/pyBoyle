#!/usr/bin/python

import os
import dill
import time

"""
Simulation Logging Framework

A custom framework to handle the creation
of log files, data files during simulation for
further analysis.
"""


class BoyleOutput(object):
    """Resultant export object from the simulation."""
    def __init__(self, exp_name, model):
        self._name = exp_name
        self._model_name = model
        self._creation_date = time.gmtime()
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
            processed=["time", "volume", "carbois", "carboin", "carbon",
                       "lipids", "lcfa", "protis", "protin", "amino",
                       "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                       "co2", "h2s", "zplus", "h2po4", "aminus", "deadcell",
                       "carb_degr", "amino_degr", "lipid_degr", "lcfa_degr",
                       "prop_degr", "butyr_degr", "valer_degr", "acet_degr",
                       "gfnh3", "gfch4", "gfco2", "gfh2s", "gasrate"],
            ph=["time", "ph"]
        )
        try:
            base_path = "./logs"
            if not os.path.exists(base_path):
                os.mkdir(base_path)
            self._path = os.path.join(base_path, self._name)
            if not os.path.exists(self._path):
                os.mkdir(os.path.join(base_path, self._name))
        except:
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
        except:
            raise

    def _update(self, attrname, value):
        """Update the attribute if the value exists"""
        try:
            getattr(self, attrname)
        except AttributeError as e:
            setattr(self, attrname, [self.__available_headers.get(attrname)])
            # -- log the error in the logging file.
        else:
            getattr(self, attrname).append(value)
