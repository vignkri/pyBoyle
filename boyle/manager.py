#!/usr/bin/env python

import numpy as np
import scipy.integrate
from export import BoyleOutput
from logger import manager_logger

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
            manager_logger.error("KeyError: Check Configuration File", e)
            print("KeyError: Check configuration file. Key is missing.")
            print("Message: ", e)
            raise
        else:
            manager_logger.info("Finished setting up model. %s" % self._meta)
            self._data_output = BoyleOutput(self._meta, self._model)

    def initialize_solver(self, iname):
        """Initialize the solver for computation"""
        try:
            self._solver = scipy.integrate.ode(self._model) \
                .set_integrator(iname, **self._frame._solver)
        except:
            raise
        finally:
            self._solver.set_initial_value(y=self._frame.Initial.get("value"),
                                           t=self._initial_time)

    def post_process(self):
        for idx in reversed(list(range(1, len(self.result)))):
            self.result[idx][29:] = (self.result[idx][29:] -
                                     self.result[idx-1][29:]) / (
                self.result[idx][0] - self.result[idx-1][0]
            ) / 1000
            self.final_result = np.hstack([self.result[idx],
                                           np.array(np.sum(self.result[idx][29:]))])
            self._data_output._update("solution", self.final_result)
        # --
        self._data_output.as_pickle()
        manager_logger.info("Post processing of data finished.")

    def __solver_start(self):
        """Start the solver"""
        manager_logger.info("Starting the solver")
        self.result = []
        # --
        while self._solver.successful() and self._solver.t < self._end_time:
            y_dot = self._solver.integrate(self._solver.t + self._step,
                                           step=True)
            # self._data_output._update("result", [self._solver.t] + list(y_dot))
            row = np.hstack([np.array([self._solver.t]), y_dot])
            self.result.append(row)
        # --
        manager_logger.info("Starting post-process of result data.")
        self.post_process()

    def function_parameters(self):
        """Pass function parameters to the simulator"""
        try:
            self._solver.set_f_params(*[self._frame, self._data_output])
        except:
            raise
        finally:
            manager_logger.info("Function parameters set in model.")

    def start(self):
        manager_logger.info("Starting the simulation.")
        self._frame.regulation()
        # -- set timing by regulation settings
        self._initial_time = 0
        self._end_time = self._frame.regulation_values["tp"][0]
        self._frame.process_data(index=0)
        self.initialize_solver(iname="vode")
        self.function_parameters()
        try:
            self.__solver_start()
        except:
            raise
        finally:
            manager_logger.info("Simulation finished successfully.")
