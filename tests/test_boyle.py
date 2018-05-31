from boyle import io
from boyle import __version__

dataset = None


def test_version():
    assert __version__ == '0.5.0'


def checkConstants(ds, raw):
    assert (ds == raw).all()


def test_initialisation():
    """Simulate settings of the testing folder."""
    import numpy as np
    # --
    data_path = "data/"
    # -- Create sample dataset
    _metadata_ = {"name": "Pytesting Data",
                  "data": data_path, "description": "test description",
                  "tags": "sample, tag, setup"}
    _settings_ = {"process": "debug",
                  "constants": "standard", "step_size": 0.5,
                  "ph": {"method": "fixed"},
                  "solver": {"method": "bdf", "order": 1,
                             "nsteps": 2, "relative": 1e-4,
                             "absolute": 1e-8}}
    # --
    _config_ = {"metadata": _metadata_, "settings": _settings_}
    # --
    global dataset
    dataset = io(_config_)
    # --
    const_one = np.loadtxt("data/Const1.constant", comments="%")
    const_two = np.loadtxt("data/Const2.constant", comments="%")
    # print(dataset.Const1.get("value"))
    checkConstants(ds=dataset.Const1.get("value"), raw=const_one)
    checkConstants(ds=dataset.Const2.get("value"), raw=const_two)
    # -- Assert values stored as solver settings data
    assert dataset._ph_method == _settings_.get("ph").get("method")
    for k in _settings_.get("solver").keys():
        if k == "relative":
            dk = "rtol"
            assert dataset._solver.get(dk) == _settings_.get("solver").get(k)
        elif k == "absolute":
            dk = "atol"
            assert dataset._solver.get(dk) == _settings_.get("solver").get(k)
        else:
            assert dataset._solver.get(k) == _settings_.get("solver").get(k)
    # -- Assert values stored as simulation configuration data
    assert dataset._simulation_config.get("step") == _settings_.get("step_size")
    assert dataset._simulation_config.get("metadata") == _metadata_.get("name")
