#!/usr/bin/env python

"""
Simulation Manager

Tools to manage the simulation started using boyle. This
enables the control of the simulation and connects the model
to the main application.
"""


import numpy as np
import scipy.integrate

from boyle.core.generic import Dataset, pHvalue
from boyle.core.load import from_localpath
from boyle.core.model.standard import Standard

# GENERIC SETTINGS
STANDARD_PH = {"method": "fixed", "value": 7.5}
STANDARD_SOLVER = {"method": "bdf", "order": 1, "nsteps": 1e6,
                   "rtol": 1e-4, "atol": 1e-8}


class Manager:
    def __init__(self, source, ph=None, solver=None,
                 step_size=0.5, dump_internals=False,
                 model="standard"):
        """Initialize manager for creating a simulation

        PARAMETERS
        ----------
        path : data
        model_name : str
        config : dict
        """
        # Get data from the local path
        if isinstance(source, Dataset):
            self._frame = source
        else:
            _data = from_localpath(source)
            self._frame = Dataset(**_data)
        # -- Get model information for setting model parameters
        if model == "standard":
            self._model = Standard
        else:
            print("Unknown model requested.")
            raise(ValueError)
        # Get pH information from the standard solver if not provided
        if not ph:
            self._ph_settings = pHvalue(STANDARD_PH.get("method"),
                                        STANDARD_PH.get("value"))
        else:
            self._ph_settings = pHvalue(ph.get("method"),
                                        ph.get("value"))
        # -- Set simulation Configuration
        if not solver:
            self._solver_setting = STANDARD_SOLVER
        else:
            self._solver_setting = solver
        # -- Get simulation configuration
        self._step = step_size
        self.integrator_name = "vode"

    def initialize_solver(self, iname):
        """Initialize the solver for computation"""
        if iname == "vode":
            self._solver = scipy.integrate.ode(self._model) \
                .set_integrator(iname, **self._solver_setting)
        elif iname == "lsoda":
            self._solver = scipy.integrate.ode(self._model) \
                .set_integrator(iname, **self._solver_setting)
        else:
            e = "ValueError: Unknown Solver provide"
            raise(e)
        # --
        self._solver.set_initial_value(y=self.initial_value,
                                       t=self._initial_time)

    def start(self, dense=False, relaxed=False):
        # Create result object to store results in
        self.result = []
        # -- loop through all available time-points to generate
        # the simulation of feeding on multiple different days.
        for idx in range(0, len(self._frame.feed_payload["tp"])):
            if idx == 0:
                self._initial_time = 0
                self._end_time = self._frame.feed_payload["tp"][idx]
            else:
                # Increase timesteps by the next timestep where there is
                # feed. Previously, the data was structured to have dT instead
                # of T of feed. Therefore there was an addition setup for
                # computing end_time. That is no longer needed.
                self._initial_time = self._end_time
                self._end_time = self._frame.feed_payload["tp"][idx]
            # --
            # -- move the io-object one index to get the required data
            self._frame.move_index_for_iteration(index=idx)
            # -- get new inoculum value from the io-object
            self.initial_value = self._frame.inoculum.get("value")
            # -- initialise the solver and the details of the solver
            self.initialize_solver(iname=self.integrator_name)
            # -- set up function parameters for a particular run_no
            _args_ = [self._frame, idx, self._ph_settings]
            self._solver.set_f_params(*_args_)
            # -- start the solver with teh current run_no
            _step = dense  # Check if there is a requirement for dense output
            _relax = relaxed  # Check if there is a req. for relaxed output
            try:
                while self._solver.successful() and \
                        self._solver.t < self._end_time:
                    try:
                        y_dot = self._solver.integrate(
                            self._solver.t + self._step, step=_step,
                            relax=_relax)
                        row = np.hstack(
                            [np.array([idx, self._solver.t]), y_dot])
                        self.result.append(row)
                    except ValueError as e:
                        error = "ValueError: The pH is diverging."
                        error += " Check the substrate and the pH computation."
                        error_payload = {"result": self.result,
                                         "internals": self._frame.debug}
                        return error_payload
            except KeyboardInterrupt as e:
                er = "KEYBOARD INTERRUPT: Stopped."
                er += " Current Iteration {}".format(self._solver.t)
            self._end_time = self._solver.t
            # The result chooses the elements from the
            # start of y_dot instead of the initial value set.
            # Forcing to use the result setup is probably not useful
            self._frame.inoculum.update({"value": y_dot})
        # --
        self._frame._update("y_hat", self.result)
        return self._frame
