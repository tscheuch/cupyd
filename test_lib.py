import platform
from pathlib import Path

import flopy
import matplotlib.pyplot as plt

from georef import CoupledModel

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
