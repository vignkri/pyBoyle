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
            self._data_output = frame

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
        for idx in reversed(list(range(1, len(self.result)))):
            self.result[idx][29:] = (self.result[idx][29:] -
                                     self.result[idx-1][29:]) / (
                self.result[idx][0] - self.result[idx-1][0]
            ) / 1000
            secondary_array = np.array(np.sum(self.result[idx][29:]))
            self.final_result = np.hstack([self.result[idx],
                                           secondary_array])
            self._data_output._update("solution", self.final_result)
        # --

    def __solver_start(self):
        """Start the solver"""
        simulationLogger.info("Starting the solver")
        self.result = []
        # --
        while self._solver.successful() and self._solver.t < self._end_time:
            y_dot = self._solver.integrate(self._solver.t + self._step,
                                           step=True)
            # self._data_output._update("result", [self._solver.t] + list(y_dot))
            row = np.hstack([np.array([self._solver.t]), y_dot])
            self.result.append(row)
        # --
        simulationLogger.info("Starting post-process of result data.")
        self._frame.Initial.update({"value": self.result[-1][1:]})

    def function_parameters(self):
        """Pass function parameters to the simulator"""
        try:
            self._solver.set_f_params(*[self._frame, self._data_output])
        except:
            raise
        finally:
            simulationLogger.info("Setting function parameters for the model")

    def start(self):
        simulationLogger.info("Starting experiment simulation.")
        # -- set timing by regulation settings
        for idx in range(0, len(self._frame.regulation_values["tp"])):
            if idx == 0:
                self._initial_time = 0
                self._end_time = self._frame.regulation_values["tp"][idx]
            else:
                self._initial_time = self._end_time
                self._end_time = self._end_time + \
                    self._frame.regulation_values["tp"][idx]
            self._frame.process_data(index=idx)
            self.initial_value = self._frame.Initial.get("value")
            self.initialize_solver(iname="vode")
            self.function_parameters()
            try:
                self.__solver_start()
            except:
                raise
            finally:
                simulationLogger.info("Simulation finished successfully.")
        # --
        self.post_process()
        self._data_output.persist()
        simulationLogger.info("Post processing of data finished.")
