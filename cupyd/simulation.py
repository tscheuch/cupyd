import os
import platform
from pathlib import Path

import flopy
import numpy
import pandas
from flopy.modflow import Modflow
from pyswmm import Nodes, Simulation, Subcatchments
from pyswmm.swmm5 import PYSWMMException

from cupyd.georef import CoupledModel





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

    @property
    def _cupled_subcatchments(self):
        return [subcatchment for subcatchment in Subcatchments(self)]

    @property
    def _cupled_storage_units(self):
        return [node for node in Nodes(self) if node.is_storage()]

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

    def before_step(self):
        """Get Before Step Callback.

        :return: Callbacks
        """
        self.set_subcatchments_cumulative_infiltration()
        self.set_storage_unites_cumulative_exfil_loss()
        # TODO: This step advance implementation con be changed by the library user and that would affect the result.
        # Also, the step advance might not always be the same!
        self.step_advance(3600 * 24)  # 1 day
        return self._callbacks["before_step"]

    def after_step(self):
        """Get After Step Callback.

        Here we override the pySWMM `Simulation` `after_step` function in order
        to call directly the coupling logic first. Any after step callback
        set is going to be called by the `_execute_coupling_logic` method.

        Check `_execute_coupling_logic` doc strings for a better understanding.

        :return: Callbacks
        """
        self._execute_coupling_logic()
        return self._callbacks["after_step"]

    def set_subcatchments_cumulative_infiltration(self):
        for subcatchment in self._cupled_subcatchments:
            self.subcatchments_cumulative_infiltration[
                subcatchment.subcatchmentid
            ] = subcatchment.statistics["infiltration"]

    def set_storage_unites_cumulative_exfil_loss(self):
        for storage_unit in self._cupled_storage_units:
            self.storage_units_cumulative_infiltration[
                storage_unit.nodeid
            ] = storage_unit.storage_statistics["exfil_loss"]

    def execute(self):
        """
        Open an input file, run SWMM, then close the file.

        Examples:

        >>> sim = PYSWMM(r'\\test.inp')
        >>> sim.execute()
        """
        for _ in self:
            pass
        self.report()
        self.close()
        # self._model.swmmExec()

    def _execute_coupling_logic(self):
        modflow_recharge_from_subcatchments = {}
        modflow_recharge_from_storage_units = {}
        for subcatchment in Subcatchments(self):
            # Delta infiltration
            modflow_recharge_from_subcatchments[subcatchment.subcatchmentid] = (
                subcatchment.statistics["infiltration"]
                - self.subcatchments_cumulative_infiltration[subcatchment.subcatchmentid]
            )

        for su in self._cupled_storage_units:
            # Delta infiltration
            # Do this to get the actual node as a `Storage`
            storage_unit = Nodes(self)[su.nodeid]

            modflow_recharge_from_storage_units[storage_unit.nodeid] = (
                storage_unit.storage_statistics["exfil_loss"]
                - self.storage_units_cumulative_infiltration[storage_unit.nodeid]
            )

        modflow_recharge_from_subcatchments_series = pandas.Series(
            modflow_recharge_from_subcatchments
        )
        modflow_recharge_from_storage_units_series = pandas.Series(
            modflow_recharge_from_storage_units
        )

        modflow_recharge_from_subcatchments_series = (
            modflow_recharge_from_subcatchments_series / self.subcatchment_area_dataframe
        )
        modflow_recharge_from_subcatchments_series.name = "subcatchment_recharge"
        modflow_recharge_from_subcatchments_series.index.name = "subcatchment"

        modflow_recharge_from_storage_units_series = (
            modflow_recharge_from_storage_units_series / self.storage_unit_area_dataframe
        )
        modflow_recharge_from_storage_units_series.name = "infiltration_storage_unit_recharge"
        modflow_recharge_from_storage_units_series.index.name = "infiltration_storage_unit"

        self.dataframe_with_recharges = pandas.merge(
            self._coupled_model.geo_dataframe,
            modflow_recharge_from_subcatchments_series,
            on="subcatchment",
            how="left",
        )
        self.dataframe_with_recharges = pandas.merge(
            self.dataframe_with_recharges,
            modflow_recharge_from_storage_units_series,
            on="infiltration_storage_unit",
            how="left",
        )

        # Aggregate cell recharges
        self.dataframe_with_recharges[
            "iteration_recharge"
        ] = self.dataframe_with_recharges.subcatchment_recharge.fillna(
            0
        ) + self.dataframe_with_recharges.infiltration_storage_unit_recharge.fillna(
            0
        )

        # Create MODFLOW inputs: RCH package (it doesn't take into account initial recharge)
        top_layer_recharge_matrix = (
            self.dataframe_with_recharges["iteration_recharge"]
            .fillna(0)
            .to_numpy()
            .reshape(self.nrows, self.ncols)
        )

        # TODO: MAKE IPAKCB GENERIC
        recharge_package = flopy.modflow.ModflowRch(
            self.modflow_model, nrchop=3, rech=top_layer_recharge_matrix, ipakcb=53
        )

        # Run MODFLOW
        # TODO: Improve performance by writing only necessary packages
        self.modflow_model.write_input()
        self.modflow_model.run_model(silent=True)

        # Read MODFLOW outputs
        # headfile, _, _ = self.modflow_model.load_results()
        fname = os.path.join(self._coupled_model.modflow_model._model_ws, "LLANQUIHUE.hds")
        headfile = flopy.utils.HeadFile(fname, model=self.modflow_model)
        heads = headfile.get_data()
        heads[heads == 1.0e30] = numpy.nan  # fix masked data
        heads[heads == -999.99] = numpy.nan

        # Strt next loop
        strt = heads[0]

        ibound = (
            self.dataframe_with_recharges["ibound"]
            .fillna(0)
            .to_numpy()
            .reshape(1, self.nrows, self.ncols)
        )

        bas = flopy.modflow.ModflowBas(
            self.modflow_model, ibound=ibound, strt=strt
        )  # use the head table of the last time step and bc

        # TODO: PREGUNTAR TERUCA
        # Profundidad a la que drena una columna de Modflow (sólo nos importa la top layer)
        DRN_burn_depth = 0.0
        # Global parameters needed to calculate drain conductance (see reference MODELMUSE DRN package pane)
        W = 5000  # model size (X)
        H = 5000  # model size (Y)
        x_resolution = W / self.ncols
        y_resolution = H / self.nrows
        DRN_L = x_resolution
        DRN_W = y_resolution
        DRN_M = 1
        DRN_K = 0.05  # m/day

        top = self.modflow_model.dis.top.array
        DTWT = (top - DRN_burn_depth) - heads[0]

        # DRN calculation
        delta_H = numpy.reshape(DTWT, len(self.dataframe_with_recharges))
        altura = numpy.reshape(heads[0], len(self.dataframe_with_recharges))
        for i in range(len(delta_H)):
            if delta_H[i] < 0:
                delta_H[i] = -delta_H[i]
            else:
                delta_H[i] = 0.0
        self.dataframe_with_recharges["Altura"] = altura
        self.dataframe_with_recharges["delta_H"] = delta_H
        self.dataframe_with_recharges["DRN_rate"] = 0.0

        # TODO: ASK TERE WHAT `DRN` IS
        self.dataframe_with_recharges["drn_cond"].fillna(0, inplace=True)
        mask = self.dataframe_with_recharges["ibound"] != -1
        self.dataframe_with_recharges.loc[mask, "DRN_rate"] = (
            self.dataframe_with_recharges["delta_H"] * self.dataframe_with_recharges["drn_cond"]
        )

        # INFLOW RATES IN SU AND JUNCTIONS

        # Inflow rates in SU:
        # node_inflow=self.dataframe_with_recharges.groupby("drn_to").sum()["DRN_rate"]
        node_inflow = self.dataframe_with_recharges.groupby("node").sum(numeric_only=True)[
            "DRN_rate"
        ]

        for node in Nodes(self):
            inflow = (
                node_inflow[node.nodeid] / 86400.0 if node.nodeid in node_inflow.index else 0
            )  # m3/s
            node.generated_inflow(inflow)

        if self._callbacks["after_step"]:
            self._callbacks["after_step"]()
