#!/usr/bin/env

"""
Input-Output Tools


Collection of tools for ingesting and exporting data
from the simulation process. The input is for creating
simulation-friendly data while the export strives to create
a dataset that is friendly to analyse and attach metadata
to.
"""


import numpy as np
from collections import namedtuple

from boyle.core.save import OUTPUT_HEADERS
from boyle.core.computations.growth import mu_max_standard
from boyle.core.computations.formula import computeHenryConstant

# order: carbis, carbin, gluc.s, prot.s, prot.in, amino, lipids,
# lcfa, hpr, hbut, hval, hac, nh4+, ch4, co2, h2s, z+, h2po4-, A-

# Standard solver settings if settings are not provided
# in configuration file
STANDARD_SOLVER_SETTINGS = {"method": "bdf", "order": 1, "nsteps": 500,
                            "relative": 1e-4, "absolute": 1e-8}

pHvalue = namedtuple("pH", "method value")


class Dataset:
    def __init__(self, **kwargs):
        """Generic Dataset Object for the simulation

        This is a generic dataset object for the simulation
        model that handles the different types of dataset
        that are passed to the simulation module. These
        objects are primarily of three different types:
            - Constants
            - Feed
            - Inoculum

        This object takes in the three different data models
        and exposes them as attributes from a single point
        for the simulation agent. The object takes keyword
        arguments which are arbitrary in order to keep the
        implementation of flexible data models.

        PARAMETERS
        ----------
        const1 : boyle.core.internals.constant.KineticConstant
            Kinetic Constant wrapper for the numpy array
            enabling the data to be accessed as a np.recarray
            for computations and indexing.
        const2 : boyle.core.internals.constant.AcidConstant
            Acid Constant wrapper for the numpy array.
        yieldc : np.array
            Data containing the yield coefficients for the
            simulation.
        feed : np.array
            The feed dataset consisting of an (n, m) matrix where
            the n is the number of days for which the simulation
            is to be run and `m` being the shape of array
            equivalent to the initial file that is defined.
        inoculum : np.array
            The inoculum file provided as a (m, 1) array where
            the m represents the required parameters for
            the simulation and the reactor composition. These
            reactor composition indices are to be preserved
            across the different arrays.

        RETURNS
        -------
        None

        REFERENCES
        ----------
        """
        self.mu_max_data = []
        self.hc_data = []
        # --
        for kw in kwargs:
            payload = {"name": kw, "value": kwargs.get(kw)}
            setattr(self, kw, payload)

        # Perform required computations for processing input
        # data. Each constant file requires a series of processing
        # and the client data requires a different series of processing
        # in order to make sure the data is usable for the simulation.
        cone_load = self.Const1.get("value").get_payload()
        self.Const1.update(dict(params=cone_load))

        # -- process client data
        self.__process_client_data()

    def _update(self, attrname, value):
        """Internal function for updating parameters

        PARAMETERS
        ----------
        attrname : str
            The name of the attribute for which the value
            has to be updated.
        value : any
            The value that needs to be updated for the
            given attribute. Currently, the attributes are
            fixed to the data loaded from the particular
            values in the dataset.

        RETURNS
        -------
        None
        """
        try:
            getattr(self, attrname)  # Try to get the attribute of given name
        except AttributeError as e:
            # If there is an attribute error, catch it
            # Check if the attribute exists in the headers file
            # if it exists, then get that value to be updated.
            if attrname in OUTPUT_HEADERS.keys():
                setattr(self, attrname, [OUTPUT_HEADERS.get(attrname)])
            else:
                # The Set the attribute with the value if not found
                # in the output_headers file.
                setattr(self, attrname, value)
        else:
            # The below is for setting list attributes with an
            # additional value.
            getattr(self, attrname).append(value)

    def __process_client_data(self):
        """Process client data to get feed and inoculum

        An internal function that processes the client data
        to create the relevant data objects for the simulation
        to utilise.  The objects that are created here are:
            a. feed information
            b. inoculum information

        The feed information is taken from the numpy array and
        converted to a dictionary called the feed payload. The
        inoculum is concatenated with the additional 4 column
        values for the gas data.
        """

        time_periods = self.feed.get("value")[:, 0]

        # TODO: Handle single day feed elegantly. Right now it
        # will break because there is no elegant solution to solve
        # 1-day feeding solution.
        if time_periods.shape[0] < 1:
            e_ = "IO: Time period provided for 1 day of feed."
            er = e_ + " Multi-day feeding is only supported not 1-day feed."
            raise ValueError(er)

        # -- Get other specific information
        temperatures = self.feed.get("value")[:, 1]
        # Divide flows by 24 to get the hourly input flow rate for
        # dataset
        flows = self.feed.get("value")[:, [2, 3]] / 24
        # Multiply the feed concentration with the actual flow
        # rate to get the flow magnitude
        substrate = flows[:, 0].reshape(-1, 1) * \
            self.feed.get("value")[:, 4:]

        # -- feed payload
        self.feed_payload = dict(tp=time_periods, temp=temperatures,
                                 flows=flows, substrates=substrate)
        # Work with Inoculum files
        _extension = np.concatenate((self.inoculum.get("value"),
                                     np.zeros(4, )))
        _value = {"value": _extension}
        self.inoculum.update(_value)

    def __recompute_mu_max(self, temp):
        """Compute Temperature Dependent Constants

        The mu-max growth values are updated every run with
        the iteration temperature provided by the feed data
        from the client.

        PARAMETERS
        ----------
        temp : float / int
            The temperature values are obtained as float
            or integer from the feed dataset from the iteration
            rows.
        """
        payload, mmax = mu_max_standard(self.Const1.get("value"),
                                        temp=temp)
        setattr(self, "mu_max", mmax)
        # --
        self.mu_max.update(payload)
        # --
        self.mu_max_data.append(payload.get("params"))

    def __recompute_hconstants(self, temp):
        """Compute henry constants

        On every iteration, the Henry Constants are recomputed using
        the iteration temperature as obtained from the feed
        dataset.

        PARAMETERS
        ----------
        temp : float / int
            The temperature values are obtained as float
            or integer from the feed dataset from the iteration
            rows.
        """
        hc, henry_c = computeHenryConstant(arr=self.Const2.get("value"),
                                           temp=temp)
        setattr(self, "henry_constants", hc)
        self.hc_data.append(henry_c)

    def move_index_for_iteration(self, index):
        """Process the imported dataset and update the values

        The function is utilised to move the index of the iteration
        and recompute the Kinetic Constant and Acid Constant
        parameters from the data sources.

        PARAMETERS
        ----------
        index : int
            The iterative index mapping to the row-index for
            the different data sources for computing the
            required parameters of input for the simulation.
        """
        # -- substrate flow information
        temp = self.feed_payload["temp"][index]
        self.flow_in = self.feed_payload["flows"][index, 0]
        self.flow_out = self.feed_payload["flows"][index, 1]
        self.substrate_flow = self.feed_payload["substrates"][index]
        # -- Update functions
        self.__recompute_mu_max(temp=temp)
        self.__recompute_hconstants(temp=temp)


class SimulationResult(object):
    def __init__(self, hdfile):
        """Simulation Result Object

        Object wrapper for the HDF5 output file to allow simplified
        access to the output and input groups and teh datasets
        stored in the groups.

        PARAMETERS
        ----------
        hdfile : HDF5 Object
            An open HDF5-object obtained from the h5py.File function
            from the python wrapper for the HDF5 file object.
        """
        self._file = hdfile

    def _access_group(self, _name_):
        """Access the data groups in the HDF5 file"""
        return self._file.get(_name_)

    def getHeaders(self, header_type):
        """Access the header strings for datasets"""
        return self._access_group("Headers").get(header_type).value

    def getDataset(self, category):
        """Access the output datasets from the simulation"""
        return self._access_group("Output").get(category).value
