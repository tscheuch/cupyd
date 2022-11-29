import os
import platform
import time
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import numpy
import pandas
from pyswmm import Links, Nodes, Subcatchments

from cupyd.georef import CoupledModel
from cupyd.simulation import CoupledSimulation

# WORKSPACES
ROOT_DIRECTORY = Path(__file__).resolve().parent
LLANQUIHUE = ROOT_DIRECTORY / "llanquihue"
MODFLOW_WORKSPACE = LLANQUIHUE / "MODFLOW"
SWMM_WORKSPACE = LLANQUIHUE / "SWMM"
SWMMGIS_DIRECTORY = SWMM_WORKSPACE / "GIS"

# MODFLOW
MODFLOW_MODEL_NAME = "LLANQUIHUE.nam"
MODFLOW_VERSION = "mfnwt"
MODFLOW_EXECUTABLE = "mfnwt.exe" if platform.system() == "Windows" else "mfnwt"

modflow_model = flopy.modflow.Modflow.load(
    MODFLOW_MODEL_NAME,
    version=MODFLOW_VERSION,
    exe_name=ROOT_DIRECTORY / MODFLOW_EXECUTABLE,
    model_ws=MODFLOW_WORKSPACE,
)

# SWMM
subcatchments = SWMMGIS_DIRECTORY / "SWMM_S.shp"
storage_units = SWMMGIS_DIRECTORY / "SWMM_SU.shp"
nodes = SWMMGIS_DIRECTORY / "SWMM_nodes_zones.shp"
coupled_model = CoupledModel(modflow_model, subcatchments, storage_units, nodes)

t1 = time.time()

with CoupledSimulation(
    coupled_model=coupled_model,
    coupled_data=None,
    inputfile=str(SWMM_WORKSPACE / "Llanquihue_base.inp"),
) as sim:

    nrows = sim.modflow_model.dis.nrow
    ncols = sim.modflow_model.dis.ncol
    # sim.step_advance(86400) # 1 day
    sim.step_advance(3600)  # 1 hour
    hours = 0
    subcatchments = [subcatchment for subcatchment in Subcatchments(sim)]
    storage_units = [node for node in Nodes(sim) if node.is_storage()]
    junctions = [node for node in Nodes(sim) if node.is_junction()]
    conduits = [link for link in Links(sim) if link.is_conduit()]
    subcatchments_cumulative_infiltration = {
        subcatchment.subcatchmentid: 0 for subcatchment in subcatchments
    }
    storage_units_cumulative_infiltration = {
        storage_unit.nodeid: 0 for storage_unit in storage_units
    }
    sim._coupled_model.geo_dataframe["area"] = sim._coupled_model.geo_dataframe.apply(
        lambda row: row.geometry.area, axis=1
    )
    subcatchment_area_dataframe = sim._coupled_model.geo_dataframe.groupby("subcatchment").sum(
        numeric_only=True
    )["area"]
    storage_unit_area_dataframe = sim._coupled_model.geo_dataframe.groupby(
        "infiltration_storage_unit"
    ).sum(numeric_only=True)["area"]
    print("STARTING SIMULATION")
    for step in sim:
        hours += 1

        if hours % 24 == 0:
            print("TIME: ", sim.current_time)

            modflow_recharge_from_subcatchments = {}
            modflow_recharge_from_storage_units = {}
            for subcatchment in Subcatchments(sim):
                # Delta infiltration
                modflow_recharge_from_subcatchments[subcatchment.subcatchmentid] = (
                    subcatchment.statistics["infiltration"]
                    - subcatchments_cumulative_infiltration[subcatchment.subcatchmentid]
                )
                subcatchments_cumulative_infiltration[
                    subcatchment.subcatchmentid
                ] = subcatchment.statistics["infiltration"]

            for su in storage_units:
                # Delta infiltration
                # Do this to get the actual node as a `Storage`
                storage_unit = Nodes(sim)[su.nodeid]

                modflow_recharge_from_storage_units[storage_unit.nodeid] = (
                    storage_unit.storage_statistics["exfil_loss"]
                    - storage_units_cumulative_infiltration[storage_unit.nodeid]
                )
                storage_units_cumulative_infiltration[
                    storage_unit.nodeid
                ] = storage_unit.storage_statistics["exfil_loss"]

            modflow_recharge_from_subcatchments_series = pandas.Series(
                modflow_recharge_from_subcatchments
            )
            modflow_recharge_from_storage_units_series = pandas.Series(
                modflow_recharge_from_storage_units
            )

            modflow_recharge_from_subcatchments_series = (
                modflow_recharge_from_subcatchments_series / subcatchment_area_dataframe
            )
            modflow_recharge_from_subcatchments_series.name = "subcatchment_recharge"
            modflow_recharge_from_subcatchments_series.index.name = "subcatchment"

            modflow_recharge_from_storage_units_series = (
                modflow_recharge_from_storage_units_series / storage_unit_area_dataframe
            )
            modflow_recharge_from_storage_units_series.name = "infiltration_storage_unit_recharge"
            modflow_recharge_from_storage_units_series.index.name = "infiltration_storage_unit"

            dataframe_with_recharges = pandas.merge(
                sim._coupled_model.geo_dataframe,
                modflow_recharge_from_subcatchments_series,
                on="subcatchment",
                how="left",
            )
            dataframe_with_recharges = pandas.merge(
                dataframe_with_recharges,
                modflow_recharge_from_storage_units_series,
                on="infiltration_storage_unit",
                how="left",
            )

            # Aggregate cell recharges
            dataframe_with_recharges[
                "iteration_recharge"
            ] = dataframe_with_recharges.subcatchment_recharge.fillna(
                0
            ) + dataframe_with_recharges.infiltration_storage_unit_recharge.fillna(
                0
            )
            # if hours % 720 == 0:
            #     dataframe_with_recharges.plot(column="iteration_recharge")
            #     print("PLOTTING RECHARGE")
            #     plt.show()

            # Create MODFLOW inputs: RCH package (it doesn't take into account initial recharge)
            top_layer_recharge_matrix = (
                dataframe_with_recharges["iteration_recharge"]
                .fillna(0)
                .to_numpy()
                .reshape(nrows, ncols)
            )

            # TODO: MAKE IPAKCB GENERIC
            recharge_package = flopy.modflow.ModflowRch(
                sim.modflow_model, nrchop=3, rech=top_layer_recharge_matrix, ipakcb=53
            )

            # Run MODFLOW
            # TODO: Improve performance by writing only necessary packages
            sim.modflow_model.write_input()
            sim.modflow_model.run_model(silent=True)

            # Read MODFLOW outputs
            # headfile, _, _ = sim.modflow_model.load_results()
            fname = os.path.join(MODFLOW_WORKSPACE, "LLANQUIHUE.hds")
            headfile = flopy.utils.HeadFile(fname, model=sim.modflow_model)
            heads = headfile.get_data()
            heads[heads == 1.0e30] = numpy.nan  # fix masked data
            heads[heads == -999.99] = numpy.nan

            # Strt next loop
            strt = heads[0]

            ibound = (
                dataframe_with_recharges["ibound"].fillna(0).to_numpy().reshape(1, nrows, ncols)
            )

            bas = flopy.modflow.ModflowBas(
                sim.modflow_model, ibound=ibound, strt=strt
            )  # use the head table of the last time step and bc

            # TODO: PREGUNTAR TERUCA
            # Profundidad a la que drena una columna de Modflow (sólo nos importa la top layer)
            DRN_burn_depth = 0.0
            # Global parameters needed to calculate drain conductance (see reference MODELMUSE DRN package pane)
            W = 5000  # model size (X)
            H = 5000  # model size (Y)
            x_resolution = W / ncols
            y_resolution = H / nrows
            DRN_L = x_resolution
            DRN_W = y_resolution
            DRN_M = 1
            DRN_K = 0.05  # m/day

            top = sim.modflow_model.dis.top.array
            DTWT = (top - DRN_burn_depth) - heads[0]

            # DRN calculation
            delta_H = numpy.reshape(DTWT, len(dataframe_with_recharges))
            altura = numpy.reshape(heads[0], len(dataframe_with_recharges))
            for i in range(len(delta_H)):
                if delta_H[i] < 0:
                    delta_H[i] = -delta_H[i]
                else:
                    delta_H[i] = 0.0
            dataframe_with_recharges["Altura"] = altura
            dataframe_with_recharges["delta_H"] = delta_H
            dataframe_with_recharges["DRN_rate"] = 0.0

            # TODO: ASK TERE WHAT `DRN` IS
            dataframe_with_recharges["drn_cond"].fillna(0, inplace=True)
            mask = dataframe_with_recharges["ibound"] != -1
            dataframe_with_recharges.loc[mask, "DRN_rate"] = (
                dataframe_with_recharges["delta_H"] * dataframe_with_recharges["drn_cond"]
            )

            # INFLOW RATES IN SU AND JUNCTIONS

            # Inflow rates in SU:
            # node_inflow=dataframe_with_recharges.groupby("drn_to").sum()["DRN_rate"]
            node_inflow = dataframe_with_recharges.groupby("node").sum(numeric_only=True)[
                "DRN_rate"
            ]
            # if hours % 720 == 0:
            #     t2 = time.time()
            #     print(t2 - t1)
            #     t1 = time.time()
            #     print("PLOTTING DRN RATE")
            #     dataframe_with_recharges.plot(column="DRN_rate")
            #     plt.show()

            for node in Nodes(sim):
                inflow = (
                    node_inflow[node.nodeid] / 86400.0 if node.nodeid in node_inflow.index else 0
                )  # m3/s
                node.generated_inflow(inflow)


gdf_final = coupled_model.geo_dataframe
print(gdf_final.head())

gdf_final.plot(column="subcatchment", legend=True)
gdf_final.plot(column="infiltration_storage_unit", legend=True)
gdf_final.plot(column="node", legend=True)
gdf_final.plot(column="drn_cond", legend=True)


def func(row):
    if row.drn_cond > 0:
        return row.node
    return


gdf_final["cell_drn"] = gdf_final.apply(func, axis=1)
# gdf_final.loc[gdf_final["drn_cond"] == np.nan, "cell_drn"] = 1
gdf_final.plot(column="cell_drn", legend=True, cmap="Spectral")
print(gdf_final.head())
# gdf_final.plot(column="SU", legend=True)
# gdf_final.subcatchment.plot()
plt.show()
