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
    assert __version__ == '0.6.1'


def checkConstants(ds, raw):
    assert (ds == raw).all()


def test_initialisation():
    """Simulate settings of the testing folder."""
    # --
    data_path = "data/"
    # -- Create sample dataset
    _metadata_ = {"name": "Pytesting Data",
                  "data": data_path, "description": "test description",
                  "tags": "sample, tag, setup"}
    _settings_ = {"process": "debug",
                  "constants": "standard", "step_size": 0.5,
                  "ph": {"method": "fixed", "value": 8},
                  "solver": {"method": "bdf", "order": 1,
                             "nsteps": 2, "relative": 1e-4,
                             "absolute": 1e-8}}
    # --
    phValue = namedtuple("pH", "method value")
    _ph_settings = phValue(_settings_.get("ph").get("method"),
                           _settings_.get("ph").get("value"))
    _config_ = {"metadata": _metadata_, "settings": _settings_}
    # --
    global dataset
    manager = Manager(_config_)
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
    for k in _settings_.get("solver").keys():
        if k == "relative":
            dk = "rtol"
            assert manager._solver_setting.get(dk) == \
                _settings_.get("solver").get(k)
        elif k == "absolute":
            dk = "atol"
            assert manager._solver_setting.get(dk) == \
                _settings_.get("solver").get(k)
        else:
            assert manager._solver_setting.get(k) == \
                _settings_.get("solver").get(k)
    # -- Assert values stored as simulation configuration data
    assert manager._step == _settings_.get("step_size")
    assert manager._meta.get("name") == _metadata_.get("name")
    # -- load reference output
    _path = "data/reference.hdf5"
    reference = SimulationResult(load.fromHDF5(_path))
    result = manager.start()
    # -- Check if the data is within thresholds
    reference_result = interpolateData(
        np.asarray(reference.getDataset("debug_solution")), time_step=1)
    test_result = interpolateData(np.asarray(result.y_hat), time_step=1)
    testing.assert_allclose(reference_result[1:], test_result[1:],
                            rtol=1e-5, atol=1e-5)


def test_loadResult():
    """Load result data and check information"""
    _path = "data/reference.hdf5"
    imported_result = SimulationResult(load.fromHDF5(_path))
    assert isinstance(imported_result.getDataset("debug"), np.ndarray)
    assert isinstance(imported_result.getDataset("solution"), np.ndarray)
    assert isinstance(imported_result.getDataset("debug_solution"), np.ndarray)
