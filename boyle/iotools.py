#!/usr/bin/env

import os
import time
import numpy as np

from boyle.save import to_hdf5, OUTPUT_HEADERS
from boyle.tools.logger import simulationLogger
from boyle.utility import load_constants, load_client_data

"""
Input-Output Tools


Collection of tools for ingesting and exporting data
from the simulation process. The input is for creating
simulation-friendly data while the export strives to create
a dataset that is friendly to analyse and attach metadata
to.
"""

# order: carbis, carbin, gluc.s, prot.s, prot.in, amino, lipids,
# lcfa, hpr, hbut, hval, hac, nh4+, ch4, co2, h2s, z+, h2po4-, A-

# Define globals for files required
SUPPORTED_EXTENSIONS = ("npy", "npz", "constant")
# Standard solver settings if settings are not provided
# in configuration file
STANDARD_SOLVER_SETTINGS = {"method": "bdf", "order": 1, "nsteps": 500,
                            "relative": 1e-4, "absolute": 1e-8}


class io:
    def __init__(self, path):
        """Set up Frame for setting up process information"""
        self.mu_max_data = []
        self.hc_data = []
        # -- Created metadata in time of simulation
        self.__creation_time = time.time()

        # -- Get path to the correct folder for accessing the files
        if not os.path.exists(path):
            _e = "Folder does not exist"
            simulationLogger.critical(_e)
            raise IOError(_e)
        else:
            folder = path

        items = list(os.walk(folder))
        fldr, lst, files = items[0]
        # -- Create tuple of file names to load and handle
        for _f in files:
            if _f.endswith(SUPPORTED_EXTENSIONS):
                name = _f.split(".")[0]
                file_path = os.path.join(fldr, _f)
                if not (name.startswith("feed") or name.startswith("inoculum")):
                    self.__set(_text=name, _path=file_path, _dType="constant")
                else:
                    self.__set(_text=name, _path=file_path, _dType="numpy")
        # The above code generates a tuple of file_name and file_path
        # for setting the correct attributes for the program.

        # The output headers and the output folders are to be
        # created so that the results can be persisted quickly
        output_folder_ = os.path.join(folder, "output")
        if not os.path.exists(output_folder_):
            os.mkdir(output_folder_)
            simulationLogger.warn("Created missing output folder")

        self._path = output_folder_

        # Perform required computations for processing input
        # data. Each constant file requires a series of processing
        # and the client data requires a different series of processing
        # in order to make sure the data is usable for the simulation.
        self.__process_constants()

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

    def __set(self, _text, _path, _dType):
        if _dType == "constant":
            if not _path.endswith("constant"):
                # Raise error if the file is found to have an extension
                # that is not valid for the constants file.
                _e = "IOERROR: File Format mismatch for {}.".format("constant")
                er_ = _e + " Given file is {}".format(_path)
                simulationLogger.error(_e)
                raise IOError(er_)
            # Create the required dictionary to be loaded at the
            # attribute name and set up the data payload
            payload = {"path": _path, "value": load_constants(_path)}
            # -- set attribute to self with the above payload
            setattr(self, _text, payload)

        elif _dType == "numpy":
            # Raise error if the file is found to have an extension
            # that is not valid for a numpy file.
            if not _path.endswith(("npy", "npz")):
                _e = "IOERROR: File Format mismatch for {}.".format("numpy")
                er_ = _e + " Given file is {}".format(_path)
                simulationLogger.error(_e)
                raise IOError(er_)
            # Create the required dictionary to be loaded at the
            # attribute name and set up the data payload
            payload = {"path": _path, "value": load_client_data(_path)}
            # -- set attribute to self with the above payload
            setattr(self, _text, payload)

    def __process_constants(self):
        const1 = self.Const1.get("value")
        const_one_payload = dict(
            kd0=0.05, ks=const1[2:, 5], ks_nh3=const1[2:, 6],
            pk_low=const1[2:, 9], pk_high=const1[2:, 10],
            ki_carbon=const1[0, 7], ki_prot=const1[1, 7],
            ki_hac_hpr=const1[6, 7], ki_hac_hbut=const1[7, 7],
            ki_nh3_hac=const1[9, 7], ki_hac_hval=const1[8, 7],
            ki_lcfa=const1[2:, 8]
        )
        self.Const1.update(dict(params=const_one_payload))

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
        const1 = self.Const1.get("value")
        mu_max = np.zeros((10, 1))
        mu_max_t0 = np.zeros((10, 1))
        # --
        for idx in range(0, 10):
            mu_max_t0[idx] = const1[idx, 0]
            alpha = const1[idx, 1]
            t0 = const1[idx, 2]
            t_opt = const1[idx, 3]
            t_max = const1[idx, 4]
            if temp < t_opt:
                mu_max[idx] = mu_max_t0[idx] + alpha * (temp - t0)
            else:
                mu_max[idx] = (mu_max_t0[idx] + alpha * (t_opt - t0)) * \
                    (t_max - temp) / (t_max - t_opt)
        # --
        setattr(self, "mu_max", dict(value=mu_max))
        # --
        payload = {"k0_carbon": mu_max[0, 0], "k0_prot": mu_max[1, 0],
                   "mu_max_t0": mu_max_t0[2:, ], "mu_max": mu_max[2:]}
        self.mu_max.update(dict(
            params=payload))
        # --
        self.mu_max_data.append(payload)
        simulationLogger.info("Process variable mu_max created.")

    def __recompute_hconstants(self, temp):
        """Compuate henry constants"""
        # -- delta tempature
        const2 = self.Const2.get("value")
        delta_temp = temp - const2[:, 1]
        # Uses the following columns: X(T0), a, b, c
        # This does the equation hc = T + dt * a + dt^2 * b + dt^3 * c
        henry_constants = const2[:, 0] + delta_temp * const2[:, 2] + \
            delta_temp**2 * const2[:, 3] + delta_temp**3 * const2[:, 4]
        # --
        hc = dict(
            k_h=henry_constants[[5, 7, 8, 11]],  # K_H results
            # -- log inverse values
            ka1_lcfa=10**(-henry_constants[0]),  # LCFA
            ka_nh4=10**(-henry_constants[6]),  # NH4+
            ka_hac=10**(-henry_constants[1]),  # HAC
            ka_hpr=10**(-henry_constants[2]),  # HPr
            ka_hbut=10**(-henry_constants[3]),  # HBut
            ka_hval=10**(-henry_constants[4]),  # HVal
            ka1_co2=10**(-henry_constants[9]),  # pKa1 CO2
            ka2_co2=10**(-henry_constants[10]),  # pKa2 CO2
            ka_h2s=10**(-henry_constants[12]),  # pKa H2S
            ka_h2po4=10**(-henry_constants[13]),  # pKa H2PO4-
            kw=10**(-henry_constants[14])
        )
        setattr(self, "henry_constants", hc)
        self.hc_data.append(hc)
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

    def save_to_file(self):
        """Store the data as a pickle."""
        # Get time of output for setting file name
        output_time = time.gmtime()
        # - get identifiers from the values.
        _year = output_time.tm_year
        _month = output_time.tm_mon
        _day = output_time.tm_mday
        _hour = output_time.tm_hour
        _minute = output_time.tm_min
        # - generate file name
        file_name = "output_Y{year}M{month}D{day}_{hour}H{m}M.hdf5".format(
            year=_year, month=_month, day=_day, hour=_hour, m=_minute
        )
        # - create path to the filename
        output_path = os.path.join(self._path, file_name)
        # -- creating file by trying
        try:
            to_hdf5(path=output_path, dataset=self)
        except OSError as e:
            simulationLogger.error("OSError: Output file exists."
                                   "Creating a new output file.")
