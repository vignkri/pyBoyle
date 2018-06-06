#!/usr/bin/env python

"""
Sampling

A collection of scripts to sample
the applicable substrate database and create
the required dataset for constructing the final
composition index.
"""


from pyDOE import lhs


def sampleLHS(data, shape, _samples):
    """Compute Latin Hypercube Samples"""
    return data.ppf(lhs(shape, samples=_samples))
