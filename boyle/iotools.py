#!/usr/bin/env

import os
import time
import h5py as h5
import numpy as np

from logger import simulationLogger

"""
Input-Output Tools


Collection of tools for ingesting and exporting data
from the simulation process. The input is for creating
simulation-friendly data while the export strives to create
a dataset that is friendly to analyse and attach metadata
to.
"""


class io:
    def __init__(self, configuration):
        """Set up Frame for setting up process information"""
        # -- Experiment information for metadata
        self.__experiment_name = configuration.get("name")
        self.__description = configuration.get("description")
        self.__experiment_tags = configuration.get("tags")
        self.__experiment_status = configuration.get("status")
        # -- Experiment settings for additional information
        self.__settings = configuration.get("settings")
        self.__ph_settings = configuration.get("ph")
        self.__solver_settings = configuration.get("solver")
        # -- Created metadata in time of simulation
        self.__creation_time = time.gmtime()
        # Simulation settings
        folder = configuration.get("data")
        try:
            assert os.path.exists(folder)
            self._folder = folder
        except AssertionError as e:
            simulationLogger.critical("Data folder not found.", e)
            raise
        else:
            items = list(os.walk(folder))
            fldr, lst, files = items[0]
            _names = [item for item in files]
            _files = [os.path.join(fldr, item) for item in files]
            self.setup_inputs(names=_names, files=_files)
        # Create Output Headers
        self.__headers = dict(
            debug=["run_no", "time", "mu_1", "mu_2", "mu_3", "mu_4", "mu_5",
                   "mu_6", "mu_7", "mu_8", "ph", "volume", "carbois",
                   "carboin", "carbon", "lipids", "lcfa", "protis",
                   "protin", "amino", "nh3", "hac", "hpr", "hbut",
                   "hval", "ch4", "co2", "h2s", "zplus", "h2po4",
                   "aminus", "deadcell", "carb_degr", "amino_degr",
                   "lipid_degr", "lcfa_degr", "prop_degr", "butyr_degr",
                   "valer_degr", "acet_degr", "gfnh3", "gfch4",
                   "gfco2", "gfh2s"],
            solution=["run_no", "time", "volume", "carbois", "carboin",
                      "carbon", "lipids", "lcfa", "protis", "protin", "amino",
                      "nh3", "hac", "hpr", "hbut", "hval", "ch4",
                      "co2", "h2s", "zplus", "h2po4", "aminus", "deadcell",
                      "carb_degr", "amino_degr", "lipid_degr", "lcfa_degr",
                      "prop_degr", "butyr_degr", "valer_degr", "acet_degr",
                      "gfnh3", "gfch4", "gfco2", "gfh2s", "gasrate"],
        )
        try:
            _base_path = self._folder
            output_fldr = os.path.join(_base_path, "output")
            if not os.path.exists(output_fldr):
                os.mkdir(output_fldr)
            self._path = output_fldr
        except AssertionError or OSError as e:
            simulationLogger.error("OSError: BoyleOutput Creation Module.")
            raise
        else:
            simulationLogger.info("Created IO toolkit for the system.")

    def __const1_parameters(self):
        """Update Const1 Parameters"""
        const1 = self.Const1.get("value")
        self.Const1.update(dict(
            params=dict(
                kd0=0.05, ks=const1[2:, 5], ks_nh3=const1[2:, 6],
                pk_low=const1[2:, 9], pk_high=const1[2:, 10],
                ki_carbon=const1[0, 7], ki_prot=const1[1, 7],
                ki_hac_hpr=const1[6, 7], ki_hac_hbut=const1[7, 7],
                ki_nh3_hac=const1[9, 8], ki_hac_hval=const1[8, 7],
                ki_lcfa=const1[2:, 8])
        ))
        simulationLogger.info("Process parameter space of Const1 updated.")

    def setup_inputs(self, names, files):
        """Setup all Simulation Model Inputs"""
        names_and_files = tuple(zip(names, files))
        for name, _file in names_and_files:
            if not (name.startswith("regulate") or name.startswith("Initial")):
                try:
                    setattr(self, name, {"path": _file,
                                         "value": np.loadtxt(_file, comments="%")})
                except:
                    print(_file)
                    raise
            else:
                setattr(self, name, {"path": _file,
                                     "value": np.load(_file)})
            if name == "Const1":
                const1 = self.Const1.get("value")
                self.Const1.update(dict(
                    params=dict(
                        kd0=0.05, ks=const1[2:, 5], ks_nh3=const1[2:, 6],
                        pk_low=const1[2:, 9], pk_high=const1[2:, 10],
                        ki_carbon=const1[0, 7], ki_prot=const1[1, 7],
                        ki_hac_hpr=const1[6, 7], ki_hac_hbut=const1[7, 7],
                        ki_nh3_hac=const1[9, 8], ki_hac_hval=const1[8, 7],
                        ki_lcfa=const1[2:, 8])
                ))
            elif name == "regulate":
                time_periods = self.regulate.get("value")[:, 0]
                if not isinstance(time_periods, np.ndarray):
                    time_periods = np.array(time_periods)
                temperatures = self.regulate.get("value")[:, 1]
                # TODO: Handle flow division by 24 elegantly
                # flows = self.regulate.get("value")[:, [2, 3]] / 24
                flows = self.regulate.get("value")[:, [2, 3]]
                substrate = flows[:, 0].reshape(-1, 1) * \
                    self.regulate.get("value")[:, 4:]
                self.regulation_values = dict(
                    tp=time_periods, temp=temperatures,
                    flows=flows, substrates=substrate)
            elif name == "Initial":
                _extend = np.concatenate((self.Initial.get("value"),
                                          np.zeros(4, )))

                _value = {"value": _extend}
                self.Initial.update(_value)
            else:
                pass
        simulationLogger.info("Input parameters created.")

    @property
    def _ph_method(self):
        """Return ph if defined in configuration file"""
        if self.__ph_settings.get("method"):
            return self.__ph_settings.get("method")
        else:
            return "standard"

    @property
    def _solver(self):
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
            step=self.__settings.get("step_size"),
            metadata=self.__experiment_name
        )
        return _simulation_config

    def __mu_max_compute(self, temp):
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

    def __compute_hconstants(self, temp):
        """Compuate henry constants"""
        # -- delta tempature
        const2 = self.Const2.get("value")
        delta_temp = temp - const2[:, 1]
        henry_constants = const2[:, 0] + delta_temp * const2[:, 2] + \
            delta_temp**2 * const2[:, 3] + delta_temp**3 * const2[:, 4]
        # --
        hc = dict(
            k_h=henry_constants[[1, 7, 8, 11]],
            # -- log inverse values
            ka1_lcfa=10**(-henry_constants[0]),
            ka_nh4=10**(-henry_constants[2]),
            ka_hac=10**(-henry_constants[3]),
            ka_hpr=10**(-henry_constants[4]),
            ka_hbut=10**(-henry_constants[5]),
            ka_hval=10**(-henry_constants[6]),
            ka1_co2=10**(-henry_constants[9]),
            ka2_co2=10**(-henry_constants[10]),
            ka_h2s=10**(-henry_constants[12]),
            ka_h2po4=10**(-henry_constants[13]),
            kw=10**(-henry_constants[14])
        )
        setattr(self, "henry_constants", hc)
        simulationLogger.info("Process variable Henry-Constants created.")

    def process_data(self, index):
        """Process the imported dataset and update the values."""
        # -- substrate flow information
        temp = self.regulation_values["temp"][index]
        self.flow_in = self.regulation_values["flows"][index, 0]
        self.flow_out = self.regulation_values["flows"][index, 1]
        self.substrate_flow = self.regulation_values["substrates"][index]
        # -- Update functions
        self.__mu_max_compute(temp=temp)
        self.__compute_hconstants(temp=temp)

    def persist(self, status="debug"):
        """Store the data as a pickle."""
        otime = time.gmtime()
        file_name = "output_Y{year}M{month}D{day}_{hour}h{min}m.hdf5".format(
            year=otime.tm_year, month=otime.tm_mon, day=otime.tm_mday,
            hour=otime.tm_hour, min=otime.tm_min)
        output_path = os.path.join(self._path, file_name)
        try:
            of = h5.File(output_path, "w")
        except OSError as e:
            simulationLogger.error("OSError: Output file exists."
                                   "Creating a new output file.")
        else:
            of.attrs["Experiment"] = np.string_(self.__experiment_name)
            of.attrs["CreationTime"] = self.__creation_time
            of.attrs["FinishTime"] = time.gmtime()
            of.attrs["Description"] = self.__description
            of.attrs["Tags"] = np.string_(self.__experiment_tags)
            # --
            headers = of.create_group("Headers")
            headers["debug"] = [np.string_(item) for item in
                                self.__headers.get("debug")]
            headers["solution"] = [np.string_(item) for item in
                                   self.__headers.get("solution")]
            # --
            out_grp = of.create_group("Output")
            out_grp["debug"] = self.debug[1:]
            out_grp["solution"] = self.solution[1:]
            # -- Deprecate this:
            if self.__experiment_status == "debug":
                out_grp["debug_solution"] = self.dbgSoln[1:]
            else:
                pass
            # --
            in_grp = of.create_group("Input")
            in_grp["regulate"] = self.regulate.get("value")
            in_grp["initial"] = self.Initial.get("value")

    def _update(self, attrname, value):
        """Update the attribute if the value exists"""
        try:
            getattr(self, attrname)
        except AttributeError as e:
            setattr(self, attrname, [self.__headers.get(attrname)])
            simulationLogger.warning("BoyleOutput: attribute '%s' not found"
                                     % (attrname))
            # -- log the error in the logging file.
        else:
            getattr(self, attrname).append(value)
