#!/usr/bin/env python

"""
Constants

Object wrappers for Constants
used during the simulation model. This
allows for computations to be handled on the object
instead of in a dataset. Allowing for flexible
computation changes instead of fixed computations.
"""

import numpy as np


class KineticConstants(np.ndarray):

    def __new__(cls, input_array, info="Kinetic Constants", kd0=0.05):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.info = info
        obj.kd0 = kd0
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.info = getattr(obj, 'info', None)

    def get_payload(self):
        payload_ = dict(
            kd0=self.kd0, ks=self[2:, 5], ks_nh3=self[2:, 6],
            pk_low=self[2:, 9], pk_high=self[2:, 10],
            ki_carbon=self[0, 7], ki_prot=self[1, 7],
            ki_hac_hpr=self[6, 7], ki_hac_hbut=self[7, 7],
            ki_nh3_hac=self[9, 7], ki_hac_hval=self[8, 7],
            ki_lcfa=self[2:, 8])
        return payload_
