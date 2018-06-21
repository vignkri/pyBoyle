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
from boyle.core.internals.constant import KineticConstant, AcidConstant

# Create set of supported extensions to make sure data
# is loaded from a specific group and nothing else
SUPPORTED_EXTENSIONS = ("npy", "npz", "constant")


def from_localpath(path):
    """Load files from local file path

    Loads the different data files from a local path
    in to the Dataset generic by passing them through
    the internal data models.

    PARAMETERS
    ----------
    path : str
        The path of the folder where the data is stored
        for the simulation model.

    RETURNS
    -------
    dict ::
        A dictionary object consisting of the different
        input files for the simulation model. The dictionary
        keys are taken from the name of the files:
            - Const1
            - Const2
            - yieldc
            - feed
            - inoculum
    """
    if not os.path.exists(path):
        _e = "Folder does not exist"
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
    """Load result file from local file path

    Loads the result HDF5 file from the local path for
    analysis of the result objects. The output object

    PARAMETERS
    ----------
    path : str ::
        The path string for the HDF5 file and the data
        from the simulation object.

    RETURNS
    -------
    out_ : h5py.File ::
        The hdf5py file object containing the input and
        the output data from the simulation model.
    """
    if not os.path.exists(path):
        _e = "Folder does not exist"
        raise IOError(_e)
    else:
        _path = path
    # --
    out_ = h5.File(_path, "r")
    return out_
