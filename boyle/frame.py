#!/usr/bin/env

import os
import numpy as np

"""
Process Frame data

The process data for the simulation is imported
and the data is cleaned up and updated to run the
simulation.
"""


class Frame:
    def __init__(self, folder):
        """Set up Frame for setting up process information"""
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
        setattr(self, "mu_max", mu_max)
        # --
        # Set up Const1 Parameters
        self.c1params = dict(
            kd0=0.05, ks=const1[2:, 5], ks_nh3=const1[2:, 6],
            pk_low=const1[2:, 9], pk_high=const1[2:, 10],
            ki_carbon=const1[0, 7], ki_prot=const1[1, 7],
            ki_hac_hpr=const1[6, 7], ki_hac_hbut=const1[7, 7],
            ki_nh3_hac=const1[9, 8], ki_lcfa=const1[2:, 8]
        )
        # -- mu_max parameters
        self.mu_max_params = dict(
            k0_carbon=mu_max[0, 0], k0_prot=mu_max[1, 0],
            mu_max_t0=mu_max_t0[2:, ], mu_max=mu_max[2:]
        )
