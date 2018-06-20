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
    def __init__(self, source, ph=None, solver=STANDARD_SOLVER,
                 step_size=0.5, model="standard"):
        """Initialize manager for creating a simulation

        The manager provides interfaces to set up the solver
        and the datasets for beginning the solver. This class
        takes a set of parameters that defines the solver
        to be used along with the settings aside from the data
        to the solver.

        PARAMETERS
        ----------
        source : str / boyle.core.generic.Dataset
            The source defines the path of the dataset or the
            source of the model dataset. If a boyle.Dataset object
            is passed to the manager, the data is immediately
            assigned. If not, the path is processed for data and
            the required data object is built.
        ph : dict
            The ph dictionary is the definition of the type of
            method to be utilised by the model for performing
            the simulation. The dictionary is to be built as
            follows:
                {"method": "fixed" | "brentq" | "newton-raphson",
                 "value": 8 | None | None | None}
            The value and the method data are mapped together. If
            the method is given as 'fixed' it is required to provide
            a value to assume fixed for the whole simulation. For
            other methods, it is not required to provide the value
            for the key "value".
        solver : dict
            The solver dictionary is the definition of the settings
            for the solver to utilise. This argument is an optional
            argument in case the solver settings need to be tweaked,
            and therefore a standard sample is provided:
                {"method": "bdf", "order": 1, "nsteps": 1e6,
                 "rtol": 1e-4, "atol": 1e-8}
            The parameters depending on the solver can be obtained
            from the Scipy.Integrate documentation [1].
        step_size : float : 0.5
            Provides the step size at which the outputs are required
            from the solver. This is also optional with the default
            being 0.5 with units in `hours`.
        model : str : "standard"
            A string to define the model function being used for the
            simulation. The list of available models:
                - standard

        REFERENCES
        ----------
        [1] Documentation for the scipy.integrate module:
            https://docs.scipy.org/doc/scipy/reference/integrate.html
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
        self._solver_setting = solver
        # -- Get simulation configuration
        self._step = step_size
        self.integrator_name = "vode"

    def initialize_solver(self, iname):
        """Initialise the solver for computation

        Initialises the solver function within the scipy.integrate
        module by wrapping the `scipy.integrate.ode` function and the
        `scipy.integrate.ode.set_integrator` method. This function helps
        in wrapping the `integrator name` and the solver settings
        together.

        PARAMETERS
        ----------
        iname : str
            A string is passed in order to set the name of the
            integrator for the solver. There are two solver which
            have been tested:
                - vode
                - lsoda
            Additional solvers are yet to be implemented.[1]

        RETURNS
        -------
        None

        REFERENCES
        ----------
        [1] Documentation for the module `scipy.integrate.ode`
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
        """
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
        """Start the solver for performing simulation

        PARAMETERS
        ----------
        dense : boolean : False
            Flag for requesting dense output from the solver
            if required for exposing the internals and verifying
            solution.
        relaxed : boolean : False
            Flag for setting the relaxation parameter within the
            integrator such that the integration would go until
            t_1 >= t and then return the results. If the dense output
            is set to True, this is not referenced. [1]

        RETURNS
        -------
        boyle.core.generic.Dataset::
            Boyle.Dataset object which contains all the data provided
            at the feed, the constants utilised and the results of the
            simulation. The results are stored under the attributes:
                - `y_hat` is the output of the scipy.integrate() function
                for each specified time point.
                - 'debug' is the internal output from the integration
                function with the datasets that are not obtained through
                the result object. This contains data for the pH and the
                BMP computation.

        REFERENCES
        ----------
        """
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
