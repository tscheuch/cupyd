import platform
import warnings
from pathlib import Path
from typing import List

warnings.simplefilter(action="ignore", category=FutureWarning)

import flopy
import matplotlib.pyplot as plt
import numpy
import pandas
from pyswmm import Links, Nodes, Simulation, Subcatchments

from georef import CoupledModel
from simulation.simulation import CoupledSimulation

MODFLOW_WORKSPACE = Path(__file__).resolve().parent
print(MODFLOW_WORKSPACE)
MODFLOW_WORKSPACE = MODFLOW_WORKSPACE / "llanquihue" / "MODFLOW"

modelname = "LLANQUIHUE.nam"
if platform.system() == "Windows":
    exe_name = "mfnwt.exe"
if platform.system() == "Darwin":
    exe_name = "mfnwt"


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
    for step in sim:
        hours += 1
        # print(sim.current_time)
        subcatchment_time_serie = sim.results.subcatchments_time_series_result
        storage_units_time_serie = sim.results.storage_units_time_series_result
        junctions_time_serie = sim.results.junctions_time_series_result
        conduits_time_serie = sim.results.conduits_time_series_result
        for subcatchment in subcatchments:
            new_row = {
                "Time": sim.current_time,
                "Precipitation (mm/h)": subcatchment.rainfall,
                "Evaporation (mm/d)": subcatchment.evaporation_loss,
                "Infiltration (mm/h)": subcatchment.infiltration_loss,
                "Runoff (m3/s)": subcatchment.runoff,
                "Runon (m3/s)": subcatchment.runon,
                "Cumulative Infiltration (m3)": subcatchment.statistics["infiltration"],
                "Cumulative Evaporation (m3)": subcatchment.statistics["evaporation"],
            }
            subcatchment_time_serie[
                subcatchment.subcatchmentid
            ] = subcatchment_time_serie[subcatchment.subcatchmentid].append(
                new_row, ignore_index=True
            )

        for storage_unit in storage_units:
            new_row = {
                "Time": sim.current_time,
                "Depth (m)": storage_unit.depth,
                "Head (m)": storage_unit.head,
                "Flooding (m3/s)": storage_unit.flooding,
                "Lateral inflow (m3/s)": storage_unit.lateral_inflow,
                "Total inflow (m3/s)": storage_unit.total_inflow,
                "Total outflow (m3/s)": storage_unit.total_outflow,
                "Volume (m3)": storage_unit.volume,
                "Losses (m3/s)": storage_unit.losses,
                "Cumulative Exfiltration Loss (m3)": storage_unit.storage_statistics[
                    "exfil_loss"
                ],
                "Cumulative Evaporation Loss (m3)": storage_unit.storage_statistics[
                    "evap_loss"
                ],
            }
            storage_units_time_serie[storage_unit.nodeid] = storage_units_time_serie[
                storage_unit.nodeid
            ].append(new_row, ignore_index=True)

        for junction in junctions:
            new_row = {
                "Time": sim.current_time,
                "Depth (m)": junction.depth,
                "Head (m)": junction.head,
                "Flooding (m3/s)": junction.flooding,
                "Lateral inflow (m3/s)": junction.lateral_inflow,
                "Total inflow (m3/s)": junction.total_inflow,
                "Total outflow (m3/s)": junction.total_outflow,
            }
            junctions_time_serie[junction.nodeid] = junctions_time_serie[
                junction.nodeid
            ].append(new_row, ignore_index=True)

        for conduit in conduits:
            new_row = {
                "Time": sim.current_time,
                "Depth (m)": conduit.depth,
                "Flow (m3/s)": conduit.flow,
            }
            conduits_time_serie[conduit.linkid] = conduits_time_serie[
                conduit.linkid
            ].append(new_row, ignore_index=True)

        # print(subcatchment_time_serie._elements)
        # s_SWMM_ts[i] =s_SWMM_ts[i].append(new_row, ignore_index=True)
        if hours % 24 == 0:
            print(sim.current_time)

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

            modflow_recharge_from_subcatchment_serie = modflow_recharge_from_subcatchment_serie / subcatchment_area_dataframe
            modflow_recharge_from_subcatchment_serie.name = "subcatchment"
            modflow_recharge_from_storage_units_serie = modflow_recharge_from_storage_units_serie / storage_unit_area_dataframe
            modflow_recharge_from_storage_units_serie.name = "infiltration_storage_unit"
            print(modflow_recharge_from_subcatchment_serie)
            print(modflow_recharge_from_storage_units_serie)

            print("FINN")
            print(sim._coupled_model.geo_dataframe)
            print(pandas.merge(sim._coupled_model.geo_dataframe, modflow_recharge_from_subcatchment_serie.to_frame(), on="subcatchment", how="right"))
            print(sim._coupled_model.geo_dataframe.iloc[0]["geometry"])
            print(sim._coupled_model.geo_dataframe.iloc[0]["geometry"].area)

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
