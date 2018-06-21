#!/usr/bin/env python

"""
Save Module

Collection of functions to output specific parts of
data in to the disk. This extracts data from the
process in different variables for storing on disk
or saving to disk for later.
"""

import os
import time
import h5py as h5
from numpy import string_


# Splitting Header Items
HEADER_START = ["run_no", "time"]

HEADER_DEBUG = ["mu_1", "mu_2", "mu_3", "mu_4", "mu_5", "mu_6",
                "mu_7", "mu_8", "pH", "raw_flow", "T_flow_in"]

HEADER_VOL = ["volume"]

HEADER_CORE = ["carb_ins", "carb_ine", "carb_sol",
               "prot_ins", "prot_ine", "amino", "lipids",
               "ac_lcfa", "ac_prop", "ac_buty", "ac_val", "ac_ace",
               "dg_nh4", "dg_ch4", "dg_co2", "dg_h2s", "io_z", "io_p",
               "io_a"]

HEADER_DEGRADERS = ["dead_cell", "degr_carb", "degr_amino", "degr_lipid",
                    "degr_lcfa", "degr_hprop", "degr_butyr", "degr_valer",
                    "degr_acet", "gf_nh3", "gf_ch4", "gf_co2", "gf_h2s"]

HEADER_END = ["gasrate"]

# Set a dictionary of headers that are to be set when saving
# outputs to file.
OUTPUT_HEADERS = dict(
    debug=HEADER_START + HEADER_DEBUG + HEADER_VOL + HEADER_CORE +
    HEADER_DEGRADERS,
    solution=HEADER_START + HEADER_VOL + HEADER_CORE + HEADER_DEGRADERS +
    HEADER_END
)


def to_hdf5(path, dataset, internals=False):
    """Save the dataset to hdf5 file"""
    _out_ = h5.File(path, "w")
    # --
    input_data_grp = _out_.create_group("Input")
    input_data_grp["feed"] = dataset.feed.get("value")
    input_data_grp["inoculum"] = dataset.inoculum.get("value")
    # --
    output_data_grp = _out_.create_group("Output")
    output_data_grp["debug"] = dataset.debug[1:]
    output_data_grp["solution"] = dataset.y_hat
    # --
    headers = _out_.create_group("Headers")
    headers["debug"] = [string_(item) for item in
                        OUTPUT_HEADERS.get("debug")]
    headers["solution"] = [string_(item) for item in
                           OUTPUT_HEADERS.get("solution")]
    # -- save functions for process computations
    if dataset.dump_internals:
        hc = _out_.create_group("henryconstants")
        hc["data"] = dataset.hc_data
        # -- mu_max dataset
        mu_max = _out_.create_group("GrowthData")
        for key in dataset.mu_max_data[0]:
            payload = []
            for row in dataset.mu_max_data:
                payload.append(row.get(key))
            # --
            mu_max[key] = payload


def to_file(_path, _dset):
    """Function to save the dataset as a file

    PARAMETERS
    ----------
    _path : str ::
        Path of the file where the dataset is to be saved
        for analysis after the simulation is complete.

    _dset : boyle.core.generic.Dataset ::
        Dataset object output from the boyle.Manager after
        the simulation is completed.

    RETURNS
    -------
    STATUS CODE : str ::
        Status code regarding the condition of operation
        after saving. The parameters are either `SUCCESS`
        or `FAIL`.
    """
    output_time = time.gmtime()
    # - get identifiers from the values.
    _year = output_time.tm_year
    _month = output_time.tm_mon
    _day = output_time.tm_mday
    _hour = output_time.tm_hour
    _minute = output_time.tm_min
    # - generate file name
    file_name = "output_Y{year}M{month}D{day}_{hour}H{m}M.hdf5".format(
        year=_year, month=_month, day=_day, hour=_hour, m=_minute
    )
    # - create path to the filename
    output_path = os.path.join(_path, file_name)
    # -- creating file by trying
    try:
        to_hdf5(path=output_path, dataset=_dset)
    except OSError as e:
        status_code = "FAIL"
    else:
        status_code = "SUCCESS"
    return status_code
