import platform
import time
from pathlib import Path

import flopy
import matplotlib.pyplot as plt

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
    hours = 0

    for step in sim:
        hours += 1

        if hours % 24 == 0:
            print("TIME: ", sim.current_time)

        if hours % 720 == 0:
            t2 = time.time()
            print(t2 - t1)
            t1 = time.time()
            sim.dataframe_with_recharges.plot(column="iteration_recharge")
            print("PLOTTING RECHARGE")
            plt.show()
            print("PLOTTING DRN RATE")
            sim.dataframe_with_recharges.plot(column="DRN_rate")
            plt.show()


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
