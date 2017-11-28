#!/usr/bin/env

import os
import dill
import time
import numpy as np

from logger import simulationLogger

"""
Process Frame data

The process data for the simulation is imported
and the data is cleaned up and updated to run the
simulation.

TODO: Refactor code for proper error catching.
"""


class Parameters:
    def __init__(self, configuration):
        """Set up Frame for setting up process information"""
        self.__experiment_name = configuration.get("name")
        self.__settings = configuration.get("settings")
        self.__solver_settings = configuration.get("solver")
        # TODO: Refactor regulate settings variable
        regulate_settings = configuration.get("regulate")

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
            self.import_files(names=_names, files=_files)

    @property
    def OutputFolder(self):
        return self._folder

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

    def import_files(self, names, files):
        """Import the files from the defined folder"""
        parsed = tuple(zip(names, files))
        for name, _file in parsed:
            setattr(self, name, {"path": _file,
                                 "value": np.loadtxt(_file, comments="%")})

    def regulation(self):
        """Set up solver regulation settings"""
        # -- Update initial values set at start
        # TODO: CHANGE THIS PROPERLY
        self.Initial.update({"value": np.concatenate(
            (self.Initial["value"], np.zeros(4, ))
        )})
        # --
        time_periods = self.regulate.get("value")[:, 0]
        if not isinstance(time_periods, np.ndarray):
            time_periods = np.array(time_periods)
        temperatures = self.regulate.get("value")[:, 1]
        flows = self.regulate.get("value")[:, [2, 3]] / 24
        substrate = flows[:,0].reshape(-1,1) * self.regulate.get("value")[:, 4:]
        self.regulation_values = dict(
            tp=time_periods, temp=temperatures,
            flows=flows, substrates=substrate
        )

    def process_data(self, index):
        """Process the imported dataset and update the values."""
        # -- substrate flow
        temp = self.regulation_values["temp"][index]
        self.flow_in = self.regulation_values["flows"][index, 0]
        self.flow_out = self.regulation_values["flows"][index, 1]
        self.substrate_flow = self.regulation_values["substrates"][index]
        # -- Compute Temperature Dependent Constants
        const1 = self.Const1.get("value")
        mu_max = np.zeros((10, 1))
        mu_max_t0 = np.zeros((10, 1))
        for idx in range(0, 10):
            mu_max_t0[idx] = const1[idx, 0]
            alpha = const1[idx, 1]
            t0 = const1[idx, 2]
            t_opt = const1[idx, 3]
            t_max = const1[idx, 4]
            #
            if temp < t_opt:
                mu_max[idx] = mu_max_t0[idx] + alpha * (temp - t0)
            else:
                mu_max[idx] = (mu_max_t0[idx] + alpha * (t_opt - t0)) * \
                    (t_max - temp) / (t_max - t_opt)
        setattr(self, "mu_max", dict(value=mu_max))
        # --
        # Set up Const1 Parameters
        self.Const1.update(dict(
            params=dict(
                kd0=0.05, ks=const1[2:, 5], ks_nh3=const1[2:, 6],
                pk_low=const1[2:, 9], pk_high=const1[2:, 10],
                ki_carbon=const1[0, 7], ki_prot=const1[1, 7],
                ki_hac_hpr=const1[6, 7], ki_hac_hbut=const1[7, 7],
                ki_nh3_hac=const1[9, 8], ki_hac_hval=const1[8, 7],
                ki_lcfa=const1[2:, 8])
        ))
        # -- mu_max parameters
        self.mu_max.update(dict(
            params=dict(
                k0_carbon=mu_max[0, 0], k0_prot=mu_max[1, 0],
                mu_max_t0=mu_max_t0[2:, ], mu_max=mu_max[2:]
            )
        ))
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
