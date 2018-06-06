#!/usr/bin/env python

"""
Composition

Tools to generate composition of the feeds
that are the input to the model and perform
quick analysis and samples of the composition based
on different techniques.
"""


import numpy as np
from scipy import stats


def createNormalDistribution(arr):
    """Create Normal Distribution for the Array"""
    try:
        assert arr.shape[1] == 2
    except AssertionError as e:
        print("Array does not have the correct shape")
        print("Construct a (n, 2) array with (mean, std.deviation)")
        raise(e)
    else:
        return stats.norm(loc=arr[:, 0], scale=arr[:, 1])
