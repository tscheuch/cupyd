from flopy.modflow import Modflow
from pyswmm import Simulation
from pyswmm.swmm5 import PYSWMMException

from cupyd.georef import CoupledModel

from .results import (
    JunctionsTimeSeriesResult,
    LinksTimeSeriesResult,
    StorageUnitsTimeSeriesResult,
    SubcatchmentsTimeSeriesResult,
)

SWMM_path = ""


class SimulationResults:
    def __init__(self, simulation) -> None:
        self.subcatchments_time_series_result = SubcatchmentsTimeSeriesResult(simulation)
        self.storage_units_time_series_result = StorageUnitsTimeSeriesResult(simulation)
        self.junctions_time_series_result = JunctionsTimeSeriesResult(simulation)
        self.conduits_time_series_result = LinksTimeSeriesResult(simulation)


def get_modflow_step_data():
    """From every MODFLOW step we need to retrieve:
    - Calculate each cell drained water, and give it back to its related Junction as lateral inflow.
        It is important to notice that this lateral inflow must be an aggregate of every Junction related cell on MODFLOW.
    - Store each cell water table depth (altura de agua) in order to initialize the MODFLOW model simulation with this as initial condition for the next step.
    """
    ...


def get_pyswmm_step_data():
    """From every PySWMM step we want to retreive:
    - Calculate the infiltration volume of the step per subcatchment (t+1) - t
    - Calculate the infiltration volume of the step per storage unit (t+1) - t

    Incorporate both of these elements as recharge for the MODFLOW model by adding it to the .rch package (before initializing the next step MODFLOW model)
    """
    ...


def make_something_with_the_data():
    """
    1. Calculate PySWMM data
    2. Add infiltration as
    """
    ...


class SmartSimulation(Simulation):
    """An extension of PySWMM simulation to allow passing callbacks to be executed in a fixed period
    of time.

    Args:
        Simulation (_type_): PySWMM simulation class
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._callbacks["add_on_x_seconds_before_step"] = None
        self._callbacks["on_x_seconds_after_step"] = None

    def on_x_seconds_before_step(self):
        """Get On Seconds X Before Step Callback.

        :return: Callbacks
        """
        return self._callbacks["on_x_seconds_before_step"]

    def add_on_x_seconds_before_step(self, callback, seconds: int):
        """
        Add callback function/method/object to execute before
        a simulation step, IF AND ONLY IF, that step simulation time is a multiple
        of the parameter X. Needs to be callable.

        :param func callback: Callable Object
        """
        if self._is_callback(callback):
            self._callbacks["on_x_seconds_before_step"] = callback

    def on_x_seconds_after_step(self):
        """Get On Seconds X After Step Callback.

        :return: Callbacks
        """
        return self._callbacks["on_x_seconds_after_step"]

    def add_on_x_seconds_after_step(self, callback, seconds: int):
        """
        Add callback function/method/object to execute before
        a simulation step, IF AND ONLY IF, that step simulation time is a multiple
        of the parameter X. Needs to be callable.

        :param func callback: Callable Object
        """
        if self._is_callback(callback):
            self._callbacks["on_x_seconds_after_step"] = callback


class CoupledSimulation(SmartSimulation):
    def __init__(self, coupled_model: CoupledModel, coupled_data=None, **kwargs):
        super().__init__(**kwargs)
        self._coupled_model = coupled_model
        self._coupled_data = coupled_data
        self.results = SimulationResults(self)

    @property
    def modflow_model(self) -> Modflow:
        """Property to access simulation related MODFLOW model easily.

        Returns:
            CoupledSimulation related MODFLOW model.
        """
        return self._coupled_model.modflow_model

    def _execute_callback(self, callback):
        """Runs the callback."""
        if callback:
            try:
                callback(self)
            except PYSWMMException:
                error_msg = "Callback Failed"
                raise PYSWMMException((error_msg))

    def _execute_coupling_logic(self):
        ...
