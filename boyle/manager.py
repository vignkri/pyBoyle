#!/usr/bin/env python

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
    def __init__(self, model, config):
        """Initialize manager for creating a simulation

        PARAMETERS
        ----------
        model : str
        config : dict
        """
        self._model = model
        try:
            self._initial_value = config.get("initial")
            self._initial_time = config.get("start_time")
            self._end_time = config.get("end_time")
            self._step = config.get("step")
            self._meta = config.get("metadata")
        except KeyError as e:
            manager_logger.error("KeyError: Check Configuration File", e)
            print("KeyError: Check configuration file. Key is missing.")
            print("Message: ", e)
            raise
        else:
            manager_logger.info("Finished setting up model. %s" % self._meta)
            self._data_output = BoyleOutput(self._meta)

    def initialize_solver(self, iname, i_params):
        """Initialize the solver for computation"""
        try:
            self._solver = scipy.integrate.ode(self._model) \
                .set_integrator(iname, **i_params)
        except:
            raise
        finally:
            self._solver.set_initial_value(y=self._initial_value,
                                           t=self._initial_time)

    def post_process(self, result):
        for idx in reversed(list(range(1, len(result)))):
            if result[idx][0] - result[idx-1][0] != 0:
                result[idx][29:] = (result[idx][29:] - result[idx-1][29:]) / (
                    result[idx][0] - result[idx-1][0]
                ) / 1000
            else:
                result[idx][29:] = (result[idx][29:] - result[idx-1][29:]) / (
                    result[idx-1][0]
                ) / 1000
            self._data_output._update("processed", result[idx])
        # --
        self._data_output.as_pickle()
        manager_logger.info("Post processing of data finished.")

    def __solver_start(self):
        """Start the solver"""
        manager_logger.info("Starting the solver")
        result = []
        while self._solver.successful() and self._solver.t < self._end_time:
            y_dot = self._solver.integrate(self._solver.t + self._step,
                                           step=True)
            self._data_output._update("result", [self._solver.t] + list(y_dot))
            result.append(y_dot)
        # --
        manager_logger.info("Starting post-process of result data.")
        self.post_process(result=result)

    def function_parameters(self, parameters):
        """Pass function parameters to the simulator"""
        try:
            parameters.append(self._data_output)
            self._solver.set_f_params(*parameters)
        except:
            raise
        finally:
            manager_logger.info("Function parameters set in model.")

    def start(self):
        manager_logger.info("Starting the simulation.")
        try:
            self.__solver_start()
        except:
            raise
        finally:
            manager_logger.info("Simulation finished successfully.")
