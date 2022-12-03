from flopy.modflow import Modflow
from pyswmm import Simulation
from pyswmm.swmm5 import PYSWMMException

from cupyd.georef import CoupledModel

SWMM_path = ""


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


class CoupledSimulation(Simulation):
    def __init__(self, coupled_model: CoupledModel, coupled_data=None, **kwargs):
        super().__init__(**kwargs)
        self._coupled_model = coupled_model
        self._coupled_data = coupled_data
        self.subcatchments_cumulative_infiltration = {}
        self.storage_units_cumulative_infiltration = {}
        self._coupled_model.geo_dataframe["area"] = self._coupled_model.geo_dataframe.apply(
            lambda row: row.geometry.area, axis=1
        )
        self.subcatchment_area_dataframe = self._coupled_model.geo_dataframe.groupby(
            "subcatchment"
        ).sum(numeric_only=True)["area"]
        self.storage_unit_area_dataframe = self._coupled_model.geo_dataframe.groupby(
            "infiltration_storage_unit"
        ).sum(numeric_only=True)["area"]
        self.nrows = self.modflow_model.dis.nrow
        self.ncols = self.modflow_model.dis.ncol
        # TODO: Delete this. Element to facilitate plotting
        self.dataframe_with_recharges = None

    @property
    def modflow_model(self) -> Modflow:
        """Property to access simulation related MODFLOW model easily.

        Returns:
            CoupledSimulation related MODFLOW model.
        """
        return self._coupled_model.modflow_model

    def _execute_callback(self, callback):
        """Runs the callback.

        This function overrides the pySWMM `Simulation` `_execute_callback` method
        in order to pass the simulation object to the callbacks. We do this so that
        we can execute the coupling logic as a "callback"; also because we believe
        it is useful to have it on the callbacks.
        
        """
        if callback:
            try:
                callback(self)
            except PYSWMMException:
                error_msg = "Callback Failed"
                raise PYSWMMException((error_msg))

    def after_step(self):
        """Get After Step Callback.

        Here we override the pySWMM `Simulation` `after_step` function in order
        to call directly the coupling logic first. Any after step callback
        set is going to be called by the `_execute_coupling_logic` method.

        Check `_execute_coupling_logic` doc strings for a better understanding.

        :return: Callbacks
        """
        return self._execute_coupling_logic(self)

    def _execute_coupling_logic(self):

        if self._callbacks["after_step"]:
            self._callbacks["after_step"]()
