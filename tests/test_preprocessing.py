from numpy.testing import assert_almost_equal
from scipy.stats import norm
from boyle import sampleLHS


def test_sampleLHS():
    """Test LatinHyperCubeSampling function"""
    data = norm(5, 15)
    sampled_data = (sampleLHS(data, shape=1, _samples=500))
    assert_almost_equal(data.mean(), sampled_data.mean(), 1e-4)
    assert_almost_equal(data.var(), sampled_data.var(), 1e-4)
