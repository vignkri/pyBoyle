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
from boyle.preprocessing.sampling import sampleLHS


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


def createSampledComposition(df, _method="lhs", sample_size=100):
    """Create Composition Sample using Custom Method"""
    # -- create rule based datasets
    normable = df.loc[df.STDEV > 0, :]  # sample when std.dev is avail.
    averaged = df.loc[np.logical_and(df["AVG"] != 0, df["STDEV"] == 0), :]
    # -- create indices
    idx_all = set(df.index)
    idx_norm = set(normable.index)
    idx_avg = set(averaged.index)
    idx_zeros = idx_all.difference(set.union(idx_norm, idx_avg))
    # -- create distributions
    if _method == "lhs":
        normal_ = createNormalDistribution(normable.values)
        out_sampler = sampleLHS(normal_, normable.shape[0],
                                sample_size)
        sampled_ = tuple(zip(idx_norm, out_sampler.T))
    else:
        pass
    # -- create ones and zeros mask
    ones_mask = np.ones((sample_size, len(idx_avg)))
    zero_mask = np.zeros((sample_size, len(idx_zeros))).T
    # -- create averaged samples
    averaged_ = tuple(zip(idx_avg, (averaged["AVG"].values * ones_mask).T))
    # -- create zero samples
    zeros_ = tuple(zip(idx_zeros, zero_mask))
    temp_out_ = [sampled_, averaged_, zeros_]
    out_ = [item for slist in temp_out_ for item in slist]
    return out_
