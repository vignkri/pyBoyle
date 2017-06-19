#!/usr/bin/env python

import scipy.integrate

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
        self._initial_value = config.get("initial")
        self._initial_time = config.get("start_time")
        self._end_time = config.get("end_time")
        self._step = config.get("step")
        self._meta = config.get("metadata")

    def initialize_solver(self, iname, i_params):
        """Initialize the solver for computation"""
        self._solver = scipy.integrate.ode(self._model) \
            .set_integrator(iname, **i_params)
        self._solver.set_initial_value(y=self._initial_value,
                                       t=self._initial_time)

    def __solver_start(self):
        """Start the solver"""
        while self._solver.successful() and self._solver.t < self._end_time:
            self._forecast = self._solver.integrate(self._solver.t + self._step,
                                                    step=True)

    def function_parameters(self, parameters):
        """Pass function parameters to the simulator"""
        self._solver.set_f_params(*parameters)

    def start(self):
        try:
            self.__solver_start()
        except:
            raise
