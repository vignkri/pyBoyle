#!/usr/bin/env python

"""
Growth Rate Computation

During each run, the growth is recomputed
using the mu_max values and the temperature
values. These function should be handled externally
from the Dataset object as these are manipulations
on the dataset and not just introspective / indexing
tools.
"""

import numpy as np


def mu_max_standard(arr, temp):
    """Standard function to compute the mu_max"""
    const1 = arr
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
    # setattr(self, "mu_max", dict(value=mu_max))
    # --
    payload = {"k0_carbon": mu_max[0, 0], "k0_prot": mu_max[1, 0],
               "mu_max_t0": mu_max_t0[2:, ], "mu_max": mu_max[2:]}
    return ({"params": payload}, {"value": mu_max})
