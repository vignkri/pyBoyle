#!/usr/bin/env

import os
import numpy as np

"""
Process Frame data

The process data for the simulation is imported
and the data is cleaned up and updated to run the
simulation.
"""


class Frame:
    def __init__(self, folder):
        """Set up Frame for setting up process information"""
        try:
            assert os.path.exists(folder)
            self._folder = folder
        except AssertionError as e:
            print("Folder not found.", e)
        else:
            items = list(os.walk(folder))
            fldr, lst, files = items[0]
            _names = [item for item in files]
            _files = [os.path.join(fldr, item) for item in files]
            self.import_files(names=_names, files=_files)

    def import_files(self, names, files):
        """Import the files from the defined folder"""
        parsed = tuple(zip(names, files))
        for name, _file in parsed:
            setattr(self, name, {"path": _file,
                                 "value": np.loadtxt(_file, comments="%")})
