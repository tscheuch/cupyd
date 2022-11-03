import os
import platform
import time
import warnings
from pathlib import Path

warnings.simplefilter(action="ignore", category=FutureWarning)

import flopy
import matplotlib.pyplot as plt
import numpy
import pandas
from pyswmm import Links, Nodes, Simulation, Subcatchments

from georef import CoupledModel
from simulation.simulation import CoupledSimulation

MODFLOW_WORKSPACE = Path(__file__).resolve().parent
MODFLOW_WORKSPACE2 = Path(__file__).resolve().parent
print(MODFLOW_WORKSPACE)
MODFLOW_WORKSPACE = MODFLOW_WORKSPACE / "llanquihue" / "MODFLOW"

modelname = "LLANQUIHUE.nam"
if platform.system() == "Windows":
    exe_name = "mfnwt.exe"
if platform.system() == "Darwin":
    exe_name = "mfnwt"


exe_name = Path.joinpath(MODFLOW_WORKSPACE2, exe_name)
print(exe_name, MODFLOW_WORKSPACE)

ml = flopy.modflow.Modflow.load(
    modelname, version="mfnwt", exe_name=exe_name, model_ws=MODFLOW_WORKSPACE
)

SWMM_WORKSPACE_GIS = "llanquihue/SWMM/GIS/"
subcatchment_shp = f"{SWMM_WORKSPACE_GIS}SWMM_S.shp"
storage_units_shp_file_path = f"{SWMM_WORKSPACE_GIS}SWMM_SU.shp"
nodes_shp_file_path = f"{SWMM_WORKSPACE_GIS}SWMM_nodes_zones.shp"

coupled_model = CoupledModel(
    ml, subcatchment_shp, storage_units_shp_file_path, nodes_shp_file_path
)

t1 = time.time()


def test_before_start(self):
    print("ON BEFORE START")
    print(self)
    print(self.current_time)
    print("OUT BEFORE START")


with CoupledSimulation(
    coupled_model=coupled_model,
    coupled_data=None,
    inputfile="llanquihue/SWMM/Llanquihue_base.inp",
) as sim:

    nrows = sim.modflow_model.dis.nrow
    ncols = sim.modflow_model.dis.ncol
    # sim.add_before_step(test_before_start)
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
    subcatchment_area_dataframe = sim._coupled_model.geo_dataframe.groupby(
        "subcatchment"
    ).sum()["area"]
    storage_unit_area_dataframe = sim._coupled_model.geo_dataframe.groupby(
        "infiltration_storage_unit"
    ).sum()["area"]
    print("STARTING SIMULATION")
    for step in sim:
        hours += 1

        if hours % 24 == 0:
            print("TIME: ", sim.current_time)

            modflow_recharge_from_subcatchment = {}
            modflow_recharge_from_storage_units = {}
            for subcatchment in Subcatchments(sim):
                # Delta infiltration
                modflow_recharge_from_subcatchment[subcatchment.subcatchmentid] = (
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

            modflow_recharge_from_subcatchment_serie = pandas.Series(
                modflow_recharge_from_subcatchment
            )
            modflow_recharge_from_storage_units_serie = pandas.Series(
                modflow_recharge_from_storage_units
            )

            modflow_recharge_from_subcatchment_serie = (
                modflow_recharge_from_subcatchment_serie / subcatchment_area_dataframe
            )
            modflow_recharge_from_subcatchment_serie.name = "subcatchment_recharge"
            modflow_recharge_from_subcatchment_serie.index.name = "subcatchment"

            modflow_recharge_from_storage_units_serie = (
                modflow_recharge_from_storage_units_serie / storage_unit_area_dataframe
            )
            modflow_recharge_from_storage_units_serie.name = (
                "infiltration_storage_unit_recharge"
            )
            modflow_recharge_from_storage_units_serie.index.name = (
                "infiltration_storage_unit"
            )

            dataframe_with_recharges = pandas.merge(
                sim._coupled_model.geo_dataframe,
                modflow_recharge_from_subcatchment_serie,
                on="subcatchment",
                how="left",
            )
            dataframe_with_recharges = pandas.merge(
                dataframe_with_recharges,
                modflow_recharge_from_storage_units_serie,
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
            if hours % 720 == 0:
                dataframe_with_recharges.plot(column="iteration_recharge")
                print("PLOTING RECHARGE ")
                plt.show()

            # Create MODFLOW inputs: RCH package (It does not take into accountancy initial recharge)

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
            # headfile, _,  _ = sim.modflow_model.load_results()
            fname = os.path.join(MODFLOW_WORKSPACE, "LLANQUIHUE.hds")
            headfile = flopy.utils.HeadFile(fname, model=sim.modflow_model)
            heads = headfile.get_data()
            heads[heads == 1.0e30] = numpy.nan  # fix masked data
            heads[heads == -999.99] = numpy.nan

            # Strt next loop
            strt = heads[0]

            ibound = (
                dataframe_with_recharges["ibound"]
                .fillna(0)
                .to_numpy()
                .reshape(1, nrows, ncols)
            )

            bas = flopy.modflow.ModflowBas(
                sim.modflow_model, ibound=ibound, strt=strt
            )  # use the head table of the last time step and bc

            # TODO: PREGUNTAR TERUCA
            # Profundidad a la que drena una columna de Modflow (Solo  nos importa la top layer)
            DRN_burn_depth = 0.0
            # Global parameters needed to calculate drain conductance (see reference MODELMUSE DRN package pane)
            W = 5000  # model size (X)
            H = 5000  # model size (Y)
            x_resolution = W / ncols
            y_resolution = H / nrows
            DRN_L = x_resolution
            DRN_W = y_resolution
            DRN_M = 1
            DRN_K = 0.05  # m/dia

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
            dataframe_with_recharges['drn_cond'].fillna(0, inplace=True)
            dataframe_with_recharges['DRN_rate'] = dataframe_with_recharges.apply(lambda r: r['delta_H'] * r['drn_cond'] if r['ibound'] != -1 else 0.0, axis=1)

            # INFLOW RATES IN SU AND JUNCTIONS

            # Inflow rates in SU:

            # node_inflow=dataframe_with_recharges.groupby("drn_to").sum()["DRN_rate"]
            node_inflow = dataframe_with_recharges.groupby("node").sum()["DRN_rate"]
            if hours % 720 == 0:
                t2 = time.time()
                print(t2 - t1)
                t1 = time.time()
                print("PLOTING DRN RATE ")
                dataframe_with_recharges.plot(column="DRN_rate")
                plt.show()

            for node in Nodes(sim):
                inflow = (
                    node_inflow[node.nodeid] / 86400.0
                    if node.nodeid in node_inflow.index
                    else 0
                )  # m3/s
                node.generated_inflow(inflow)


gdf_final = coupled_model.geo_dataframe

# gdf_final = georeference_models(ml, 'SWMM_inputs/shapes/SWMM_S.shp')
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
