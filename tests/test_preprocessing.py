import numpy as np
from numpy.testing import assert_almost_equal
from scipy.stats import norm
from boyle import sampleLHS
from boyle import createNormalDistribution


def test_sampleLHS():
    """Test LatinHyperCubeSampling function"""
    data = norm(5, 15)
    sampled_data = (sampleLHS(data, shape=1, _samples=500))
    assert_almost_equal(data.mean(), sampled_data.mean(), 1e-4)
    assert_almost_equal(data.var(), sampled_data.var(), 1e-4)


def test_normdistcreation():
    """Test creation of normal distribution"""
    test_array = np.vstack([np.arange(5, 10), np.arange(12, 17)]).transpose()
    assert test_array.shape[1] == 2
    # -- get mean column and variance column
    means = test_array[:, 0]
    variances = test_array[:, 1]
    # -- compute the norm values
    x = (createNormalDistribution(test_array))
    assert_almost_equal(means, x.mean(), 1e-4)
    assert_almost_equal(variances**2, x.var(), 1e-4)
