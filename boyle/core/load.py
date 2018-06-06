#!/usr/bin/env python

"""
Load

A collection of functions to load and
set up the Dataset object for performing
simulations.
"""

import os
import h5py as h5
from boyle.tools.utility import load_data
from boyle.tools.logger import simulationLogger
from boyle.core.internals.constant import KineticConstant, AcidConstant

# Create set of supported extensions to make sure data
# is loaded from a specific group and nothing else
SUPPORTED_EXTENSIONS = ("npy", "npz", "constant")


def from_localpath(path):
    """Load files from local file path"""
    if not os.path.exists(path):
        _e = "Folder does not exist"
        simulationLogger.critical(_e)
        raise IOError(_e)
    else:
        folder = path
    # Walk the folder to get the list of files
    items = list(os.walk(folder))
    fldr, lst, files = items[0]
    # -- Create tuple of file names to load and handle
    _import_data = {}
    for _f in files:
        if _f.endswith(SUPPORTED_EXTENSIONS):
            name = _f.split(".")[0]
            file_path = os.path.join(fldr, _f)
            # --
            if name == "Const1":
                value = KineticConstant(load_data(file_path))
            elif name == "Const2":
                data = load_data(file_path)
                value = AcidConstant(data)
            else:
                value = load_data(file_path)
            # -- update the importing dataset
            _import_data.update({"{}".format(name): value})
    # --
    return _import_data


def fromHDF5(path):
    """Load result file from local file path"""
    if not os.path.exists(path):
        _e = "Folder does not exist"
        simulationLogger.critical(_e)
        raise IOError(_e)
    else:
        _path = path
    # --
    out_ = h5.File(_path, "r")
    return out_
