#!/usr/bin/env python

import numpy as np
import scipy.integrate
from logger import simulationLogger

"""
Simulation Manager

Tools to manage the simulation started using boyle. This
enables the control of the simulation and connects the model
to the main application.
"""


class Manager:
    def __init__(self, model, frame):
        """Initialize manager for creating a simulation

        PARAMETERS
        ----------
        model : str
        config : dict
        """
        self._model = model
        self._frame = frame
        try:
            self._step = frame._simulation_config.get("step")
            self._meta = frame._simulation_config.get("metadata")
        except KeyError as e:
            simulationLogger.error("KeyError: Check Configuration File"
                                   "Key is missing.", e)
            raise
        else:
            simulationLogger.info("Set up experiment: '%s'" % self._meta)

    def initialize_solver(self, iname):
        """Initialize the solver for computation"""
        try:
            self._solver = scipy.integrate.ode(self._model) \
                .set_integrator(iname, **self._frame._solver)
        except:
            raise
        finally:
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

    def __solver_start(self, run_no):
        """Start the solver"""
        simulationLogger.info("Starting the solver")
        # --
        try:
            while self._solver.successful() and self._solver.t < self._end_time:
                y_dot = self._solver.integrate(self._solver.t + self._step,
                                               step=True)
                # self._data_output._update("result",
                # [self._solver.t] + list(y_dot))
                row = np.hstack([np.array([run_no, self._solver.t]), y_dot])
                self.result.append(row)
        except:
            print("Current Iteration {}".format(self._solver.t))
            raise
        # --
        # The result chooses the elements from the start of y_dot instead
        # of the initial value set. Forcing to use the result setup is probably
        # not useful
        # self._frame.Initial.update({"value": self.result[-1][2:]})
        self._frame.Initial.update({"value": y_dot})

    def function_parameters(self, run_no):
        """Pass function parameters to the simulator"""
        try:
            self._solver.set_f_params(*[self._frame, self._frame,
                                        run_no, self._frame._ph_method])
        except:
            raise
        finally:
            simulationLogger.info("Setting function parameters for the model")

    def start(self):
        simulationLogger.info("Starting experiment simulation.")
        # -- set timing by regulation settings
        self.result = []
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
            self._frame.process_data(index=idx)
            self.initial_value = self._frame.Initial.get("value")
            self.initialize_solver(iname="vode")
            self.function_parameters(run_no=idx)
            try:
                self.__solver_start(run_no=idx)
            except:
                raise
            finally:
                simulationLogger.info("Simulation finished successfully.")
        # --
        simulationLogger.info("Starting post-processing")
        # --
        self.post_process()
        self._frame.persist()
        simulationLogger.info("Post processing of data finished.")
        simulationLogger.info("Finishing up simulation. Closing.")
