#!/usr/bin/env

"""
Input-Output Tools


Collection of tools for ingesting and exporting data
from the simulation process. The input is for creating
simulation-friendly data while the export strives to create
a dataset that is friendly to analyse and attach metadata
to.
"""


import os
import time
import numpy as np

from boyle.tools.logger import simulationLogger
from boyle.core.save import to_hdf5, OUTPUT_HEADERS
from boyle.core.computations.growth import mu_max_standard
from boyle.core.computations.formula import computeHenryConstant

# order: carbis, carbin, gluc.s, prot.s, prot.in, amino, lipids,
# lcfa, hpr, hbut, hval, hac, nh4+, ch4, co2, h2s, z+, h2po4-, A-

# Standard solver settings if settings are not provided
# in configuration file
STANDARD_SOLVER_SETTINGS = {"method": "bdf", "order": 1, "nsteps": 500,
                            "relative": 1e-4, "absolute": 1e-8}


class Dataset:
    def __init__(self, **kwargs):
        """Set up Frame for setting up process information"""
        self.mu_max_data = []
        self.hc_data = []
        # --
        for kw in kwargs:
            payload = {"name": kw, "value": kwargs.get(kw)}
            setattr(self, kw, payload)
        # --
        # -- Created metadata in time of simulation
        self.__creation_time = time.time()

        # The output headers and the output folders are to be
        # created so that the results can be persisted quickly
        # output_folder_ = os.path.join(folder, "output")
        # if not os.path.exists(output_folder_):
        #     os.mkdir(output_folder_)
        #     simulationLogger.warn("Created missing output folder")

        # self._path = output_folder_

        # Perform required computations for processing input
        # data. Each constant file requires a series of processing
        # and the client data requires a different series of processing
        # in order to make sure the data is usable for the simulation.
        cone_load = self.Const1.get("value").get_payload()
        self.Const1.update(dict(params=cone_load))

        # -- process client data
        self.__process_client_data()

    def _update(self, attrname, value):
        """Update the attribute if the value exists"""
        try:
            getattr(self, attrname)  # Try to get the attribute of given name
        except AttributeError as e:
            # If there is an attribute error, catch it
            simulationLogger.warning("Attribute {} missing.".format(attrname))
            simulationLogger.info("Creating attribute {}".format(attrname))
            # Check if the attribute exists in the headers file
            # if it exists, then get that value to be updated.
            if attrname in OUTPUT_HEADERS.keys():
                setattr(self, attrname, [OUTPUT_HEADERS.get(attrname)])
            else:
                # The Set the attribute with the value if not found
                # in the output_headers file.
                setattr(self, attrname, value)
        else:
            # The below is for setting list attributes with an
            # additional value.
            getattr(self, attrname).append(value)

    def __process_client_data(self):
        """Process client data to get feed and inoculum"""

        time_periods = self.feed.get("value")[:, 0]

        # TODO: Handle single day feed elegantly. Right now it
        # will break because there is no elegant solution to solve
        # 1-day feeding solution.
        if time_periods.shape[0] < 1:
            e_ = "IO: Time period provided for 1 day of feed."
            er = e_ + " Multi-day feeding is only supported not 1-day feed."
            simulationLogger.Error(e_)
            raise ValueError(er)

        # -- Get other specific information
        temperatures = self.feed.get("value")[:, 1]
        # Divide flows by 24 to get the hourly input flow rate for
        # dataset
        flows = self.feed.get("value")[:, [2, 3]] / 24
        # Multiply the feed concentration with the actual flow
        # rate to get the flow magnitude
        substrate = flows[:, 0].reshape(-1, 1) * \
            self.feed.get("value")[:, 4:]

        # -- feed payload
        feed_payload = dict(tp=time_periods, temp=temperatures,
                            flows=flows, substrates=substrate)

        # TODO: Rename regulation values to Feed_Payload to reflect
        # real world phenomena
        self.regulation_values = feed_payload

        # Work with Inoculum files
        _extension = np.concatenate((self.inoculum.get("value"),
                                     np.zeros(4, )))
        _value = {"value": _extension}
        self.inoculum.update(_value)

    def __recompute_mu_max(self, temp):
        """Compute Temperature Dependent Constants"""
        # TODO: Clean up this function. This function should
        # return arrays instead of the currently returned
        # dictionaries. This would cause downstream issues since
        # dictionaries are complicated to be brief enough for
        # future changes.
        payload, mmax = mu_max_standard(self.Const1.get("value"),
                                        temp=temp)
        setattr(self, "mu_max", mmax)
        # --
        self.mu_max.update(payload)
        # --
        self.mu_max_data.append(payload.get("params"))
        simulationLogger.info("Process variable mu_max created.")

    def __recompute_hconstants(self, temp):
        """Compute henry constants"""
        hc, henry_c = computeHenryConstant(arr=self.Const2.get("value"),
                                           temp=temp)
        setattr(self, "henry_constants", hc)
        self.hc_data.append(henry_c)
        simulationLogger.info("Process variable Henry-Constants created.")

    def move_index_for_iteration(self, index):
        """Process the imported dataset and update the values."""
        # -- substrate flow information
        temp = self.regulation_values["temp"][index]
        self.flow_in = self.regulation_values["flows"][index, 0]
        self.flow_out = self.regulation_values["flows"][index, 1]
        self.substrate_flow = self.regulation_values["substrates"][index]
        # -- Update functions
        self.__recompute_mu_max(temp=temp)
        self.__recompute_hconstants(temp=temp)
