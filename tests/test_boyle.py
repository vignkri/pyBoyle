import numpy as np
from numpy import testing
from boyle.manager import Manager
from boyle.tools.utility import load_data
from boyle import __version__
from boyle.tools.analysis import interpolateData
from boyle import load, SimulationResult, Dataset
from collections import namedtuple

dataset = None


def test_version():
    assert __version__ == '0.6.2'


def checkConstants(ds, raw):
    assert (ds == raw).all()


def test_initialisation():
    """Simulate settings of the testing folder."""
    # --
    data_path = "data/"
    # -- Create sample dataset
    standard_ph = {"method": "fixed", "value": 7.5}
    standard_solver = {"method": "bdf", "order": 1,
                       "nsteps": 1e6, "rtol": 1e-4,
                       "atol": 1e-8}
    standard_step_size = 0.5
    # --
    phValue = namedtuple("pH", "method value")
    _ph_settings = phValue(standard_ph.get("method"), standard_ph.get("value"))
    # --
    global dataset
    manager = Manager(path=data_path)
    dataset = manager._frame
    # --
    const_one = load_data("data/Const1.constant")
    const_two = load_data("data/Const2.npy")
    # print(dataset.Const1.get("value"))
    checkConstants(ds=dataset.Const1.get("value"), raw=const_one)
    checkConstants(ds=dataset.Const2.get("value"), raw=const_two)
    # -- Assert values stored as solver settings data
    assert manager._ph_settings.method == _ph_settings.method
    assert manager._ph_settings.value == _ph_settings.value
    for k in manager._solver_setting.keys():
        assert manager._solver_setting.get(k) == \
                standard_solver.get(k)
    # -- Assert values stored as simulation configuration data
    assert manager._step == standard_step_size
    # -- load reference output
    _path = "data/reference.hdf5"
    reference = SimulationResult(load.fromHDF5(_path))
    result = manager.start()
    # -- Check if the data is within thresholds
    reference_result = interpolateData(
        np.asarray(reference.getDataset("debug_solution")),
        time_col=1, time_step=1)
    test_result = interpolateData(np.asarray(result.y_hat), time_col=1,
                                  time_step=1)
    testing.assert_allclose(reference_result[1:], test_result[1:],
                            rtol=1e-5, atol=1e-5)


def test_loadResult():
    """Load result data and check information"""
    _path = "data/reference.hdf5"
    imported_result = SimulationResult(load.fromHDF5(_path))
    assert isinstance(imported_result.getDataset("debug"), np.ndarray)
    assert isinstance(imported_result.getDataset("solution"), np.ndarray)
    assert isinstance(imported_result.getDataset("debug_solution"), np.ndarray)
