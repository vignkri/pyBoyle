#!/usr/bin/env python

"""
Simulation Manager

Tools to manage the simulation started using boyle. This
enables the control of the simulation and connects the model
to the main application.
"""


import numpy as np
import scipy.integrate
from boyle import model
from boyle.iotools import io
from boyle.tools.logger import simulationLogger


class Manager:
    def __init__(self, config, model=model.standard):
        """Initialize manager for creating a simulation

        PARAMETERS
        ----------
        model : str
        config : dict
        """
        self._model = model
        # --
        self._config = config
        self._frame = io(self._config.get("metadata").get("data"))
        # -- Get simulation configuration
        self._step = self._config.get("settings").get("step_size")
        self._meta = self._config.get("metadata")
        self._sttn = self._config.get("settings")
        # -- solver settings
        solver = self._config.get("settings").get("solver")
        self._solver_setting = dict(
            method=solver.get("method"),
            order=solver.get("order"),
            rtol=solver.get("relative"),
            atol=solver.get("absolute"),
            nsteps=solver.get("nsteps")
        )
        self._ph_settings = self._sttn.get("ph").get("method")
        simulationLogger.info("Set up experiment: '%s'" % self._meta)

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
        self._solver.set_initial_value(y=self.initial_value,
                                       t=self._initial_time)

    def post_process(self):
        """Computes the change in values

        The dy/dt is computed as the y is cumulative from the
        results. Cumulative `y` is not useful for visualisation
        as the required output is change with respect to the
        previous step.

            dy/dt = y[n] - y[n-1] / t[n] - t[n-1]
        """
        # TODO: This should be re-engineered. Data should not be
        # reversed but only the values should be subtracted with previous
        # values.
        self._frame._update("debug_solution", self.result)
        for idx in reversed(list(range(1, len(self.result)))):
            self.result[idx][29:] = (self.result[idx][29:] -
                                     self.result[idx-1][29:]) / (
                self.result[idx][1] - self.result[idx-1][1]
            ) / 1000
            # The above result[idx][1] points to the time position of
            # the index. Previously it was pointing to 0 because the
            # run_no value was not included.
            secondary_array = np.array(np.sum(self.result[idx][29:]))
            self.final_result = np.hstack([self.result[idx],
                                           secondary_array])
            self._frame._update("solution", self.final_result)
        # --

    def start(self):
        simulationLogger.info("Starting experiment.")
        # Create result object to store results in
        self.result = []
        # -- loop through all available time-points to generate
        # the simulation of feeding on multiple different days.
        for idx in range(0, len(self._frame.regulation_values["tp"])):
            simulationLogger.info("Running Simulation Frame: {}".format(idx))
            if idx == 0:
                self._initial_time = 0
                self._end_time = self._frame.regulation_values["tp"][idx]
            else:
                # Increase timesteps by the next timestep where there is
                # feed. Previously, the data was structured to have dT instead
                # of T of feed. Therefore there was an addition setup for
                # computing end_time. That is no longer needed.
                self._initial_time = self._end_time
                self._end_time = self._frame.regulation_values["tp"][idx]
            # --
            # -- move the io-object one index to get the required data
            self._frame.move_index_for_iteration(index=idx)
            # -- get new inoculum value from the io-object
            self.initial_value = self._frame.inoculum.get("value")
            # -- initialise the solver and the details of the solver
            self.initialize_solver(iname="vode")
            # -- set up function parameters for a particular run_no
            _args_ = [self._frame, idx, self._ph_settings]
            self._solver.set_f_params(*_args_)
            simulationLogger.info("Setting function parameters for the model")
            # -- start the solver with teh current run_no
            _step = True
            _relax = True
            try:
                while self._solver.successful() and self._solver.t < self._end_time:
                    y_dot = self._solver.integrate(self._solver.t + self._step,
                                                   step=_step, relax=_relax)
                    # self._data_output._update("result",
                    # [self._solver.t] + list(y_dot))
                    row = np.hstack([np.array([idx, self._solver.t]), y_dot])
                    self.result.append(row)
            except KeyboardInterrupt as e:
                er = "KEYBOARD INTERRUPT: Stopped."
                er += " Current Iteration {}".format(self._solver.t)
                simulationLogger.info(er)
                print(er)
                raise(e)
            self._end_time = self._solver.t
            # The result chooses the elements from the start of y_dot instead
            # of the initial value set. Forcing to use the result setup is probably
            # not useful
            # self._frame.inoculum.update({"value": self.result[-1][2:]})
            self._frame.inoculum.update({"value": y_dot})
            # -- Log that the simulation ended correctly.
            simulationLogger.info("Simulation finished successfully.")
        # --
        simulationLogger.info("Starting post-processing")
        # --
        self.post_process()
        # self._frame.save_to_file()
        simulationLogger.info("Post processing of data finished.")
        simulationLogger.info("Finishing up simulation. Closing.")
        return self._frame
