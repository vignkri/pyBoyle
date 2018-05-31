#!/usr/bin/env

import os
import time
import h5py as h5
import numpy as np

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

# Splitting Header Items
HEADER_START = ["run_no", "time"]

HEADER_DEBUG = ["mu_1", "mu_2", "mu_3", "mu_4", "mu_5", "mu_6",
                "mu_7", "mu_8", "pH", "raw_flow", "T_flow_in"]

HEADER_VOL = ["volume"]

HEADER_CORE = ["carb_ins", "carb_ine", "carb_sol",
               "prot_ins", "prot_ine", "amino", "lipids",
               "ac_lcfa", "ac_prop", "ac_buty", "ac_val", "ac_ace",
               "dg_nh4", "dg_ch4", "dg_co2", "dg_h2s", "io_z", "io_p",
               "io_a"]

HEADER_DEGRADERS = ["dead_cell", "degr_carb", "degr_amino", "degr_lipid",
                    "degr_lcfa", "degr_hprop", "degr_butyr", "degr_valer",
                    "degr_acet", "gf_nh3", "gf_ch4", "gf_co2", "gf_h2s"]

HEADER_END = ["gasrate"]

# Set a dictionary of headers that are to be set when saving
# outputs to file.
OUTPUT_HEADERS = dict(
    debug=HEADER_START + HEADER_DEBUG + HEADER_VOL + HEADER_CORE +
    HEADER_DEGRADERS,
    solution=HEADER_START + HEADER_VOL + HEADER_CORE + HEADER_DEGRADERS +
    HEADER_END
)


class io:
    def __init__(self, configuration):
        """Set up Frame for setting up process information"""
        # -- Experiment information for metadata
        metadata = configuration.get("metadata")
        self.__experiment_name = metadata.get("name")  # experiment name
        self.__description = metadata.get("description")  # description
        self.__experiment_tags = metadata.get("tags")  # tags for searching

        # -- Experiment settings for additional information
        settings = configuration.get("settings")
        self.__experiment_status = settings.get("process")
        self.__step_size = settings.get("step_size")
        self.__ph_settings = settings.get("ph")
        self.__solver_settings = settings.get("solver")

        # -- default to standard if solver settings is given as none
        if not self.__solver_settings:
            self.__solver_settings = STANDARD_SOLVER_SETTINGS

        # -- Created metadata in time of simulation
        self.__creation_time = time.time()

        # -- Get path to the correct folder for accessing the files
        if not os.path.exists(metadata.get("data")):
            _e = "Folder does not exist"
            simulationLogger.critical(_e)
            raise IOError(_e)

        items = list(os.walk(folder))
        fldr, lst, files = items[0]
        # -- Create tuple of file names to load and handle
        _input_data = [(_f.split(".")[0], os.path.join(fldr, _f))
                       for _f in files if _f.endswith(SUPPORTED_EXTENSIONS)]
        # The above code generates a tuple of file_name and file_path
        # for setting the correct attributes for the program.

        # The output headers and the output folders are to be
        # created so that the results can be persisted quickly
        output_folder_ = os.path.join(folder, "output")
        if not os.path.exists(output_folder_):
            os.mkdir(output_folder_)
            simulationLogger.warn("Created missing output folder")

        self._path = output_folder_

        # Start the processing of the imported data paths now
        self.setup_inputs(dataLocn=_input_data)

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

    @property
    def _ph_method(self):
        """Return pH method if defined in configuration file"""
        if self.__ph_settings.get("method"):
            return self.__ph_settings.get("method")
        else:
            # Returns standard pH method if there is no
            # defined pH method.
            return "newton-raphson"

    @property
    def _solver(self):
        """Return Solver settings for starting the simulation"""
        _solver_params = dict(
            method=self.__solver_settings.get("method"),
            order=self.__solver_settings.get("order"),
            rtol=self.__solver_settings.get("relative"),
            atol=self.__solver_settings.get("absolute"),
            nsteps=self.__solver_settings.get("nsteps")
        )
        return _solver_params

    @property
    def _simulation_config(self):
        """Simulation configuration dictionary"""
        _simulation_config = dict(
            step=self.__step_size,
            metadata=self.__experiment_name
        )
        return _simulation_config

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

    def setup_inputs(self, dataLocn):
        """Setup all Simulation Model Inputs"""
        for name, _file in dataLocn:
            if not (name.startswith("feed") or name.startswith("inoculum")):
                self.__set(_text=name, _path=_file, _dType="constant")
            else:
                self.__set(_text=name, _path=_file, _dType="numpy")
        # --

        # Perform required computations for processing input
        # data. Each constant file requires a series of processing
        # and the client data requires a different series of processing
        # in order to make sure the data is usable for the simulation.
        self.__process_constants()

        # -- process client data
        self.__process_client_data()

        # -- update logs
        simulationLogger.info("Input parameters created.")

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
        self.mu_max.update(dict(
            params=dict(
                k0_carbon=mu_max[0, 0], k0_prot=mu_max[1, 0],
                mu_max_t0=mu_max_t0[2:, ], mu_max=mu_max[2:]
            )
        ))
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
            out_ = h5.File(output_path, "w")
        except OSError as e:
            simulationLogger.error("OSError: Output file exists."
                                   "Creating a new output file.")
        else:
            out_.attrs["Experiment"] = np.string_(self.__experiment_name)
            out_.attrs["CreationTime"] = self.__creation_time
            ending_time = time.time()
            out_.attrs["FinishTime"] = ending_time
            out_.attrs["Duration"] = ending_time - self.__creation_time
            out_.attrs["Description"] = self.__description
            out_.attrs["Tags"] = np.string_(self.__experiment_tags)
            # --
            headers = out_.create_group("Headers")
            headers["debug"] = [np.string_(item) for item in
                                OUTPUT_HEADERS.get("debug")]
            headers["solution"] = [np.string_(item) for item in
                                   OUTPUT_HEADERS.get("solution")]
            # --
            output_group = out_.create_group("Output")
            output_group["debug"] = self.debug[1:]
            output_group["solution"] = self.solution[1:]
            # -- Deprecate this:
            if self.__experiment_status == "debug":
                output_group["debug_solution"] = np.array(self.debug_solution)
            else:
                pass
            # --
            input_group = out_.create_group("Input")
            input_group["regulate"] = self.feed.get("value")
            input_group["initial"] = self.inoculum.get("value")
