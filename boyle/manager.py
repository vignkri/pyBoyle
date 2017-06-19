#!/usr/bin/env python

import scipy.integrate
from export import Export
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
            self._data_exporter = Export(self._meta)

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

    def __solver_start(self):
        """Start the solver"""
        manager_logger.info("Starting the solver")
        while self._solver.successful() and self._solver.t < self._end_time:
            self._forecast = self._solver.integrate(self._solver.t + self._step,
                                                    step=True)

    def function_parameters(self, parameters):
        """Pass function parameters to the simulator"""
        try:
            parameters.append(self._data_exporter)
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
