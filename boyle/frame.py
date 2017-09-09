#!/usr/bin/env

import os
import numpy as np

"""
Process Frame data

The process data for the simulation is imported
and the data is cleaned up and updated to run the
simulation.
"""


class Parameters:
    def __init__(self, configuration):
        """Set up Frame for setting up process information"""
        settings = configuration.get("settings")
        self.__solver_settings = configuration.get("solver")
        regulate_settings = configuration.get("regulate")

        # Simulation settings
        folder = configuration.get("data")
        experiment_name = configuration.get("name")
        start_time = settings.get("t_initial")
        end_time = settings.get("t_final")
        step_size = settings.get("step_size")

        try:
            assert os.path.exists(folder)
            self._folder = folder
        except AssertionError as e:
            print("Folder not found.", e)
        else:
            items = list(os.walk(folder))
            fldr, lst, files = items[0]
            _names = [item for item in files]
            _files = [os.path.join(fldr, item) for item in files]
            self.import_files(names=_names, files=_files)

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

    def import_files(self, names, files):
        """Import the files from the defined folder"""
        parsed = tuple(zip(names, files))
        for name, _file in parsed:
            setattr(self, name, {"path": _file,
                                 "value": np.loadtxt(_file, comments="%")})

    def process_data(self, temp, flow_in, flow_out):
        """Process the imported dataset and update the values."""
        # -- Update initial values set
        self.Initial.update({"value": np.concatenate(
            (self.Initial["value"], np.zeros(4, ))
        )})
        # -- substrate flow
        self.flow_in = flow_in
        self.flow_out = flow_out
        self.substrate_flow = flow_in * self.regulate.get("value")[4:]
        # -- Compute Temperature Dependent Constants
        const1 = self.Const1.get("value")
        mu_max = np.zeros((10, 1))
        mu_max_t0 = np.zeros((10, 1))
        for index in range(0, 10):
            mu_max_t0[index] = const1[index, 0]
            alpha = const1[index, 1]
            t0 = const1[index, 2]
            t_opt = const1[index, 3]
            t_max = const1[index, 4]
            #
            if temp < t_opt:
                mu_max[index] = mu_max_t0[index] + alpha * (temp - t0)
            else:
                mu_max[index] = (mu_max_t0[index] + alpha * (t_opt - t0)) * \
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
