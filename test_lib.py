
import matplotlib.pyplot as plt
import platform

import flopy

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from georef import CoupledModel



from pathlib import Path
MODFLOW_WORKSPACE = Path(__file__).resolve().parent
print(MODFLOW_WORKSPACE)
MODFLOW_WORKSPACE = MODFLOW_WORKSPACE / 'llanquihue' / 'MODFLOW'

modelname = 'LLANQUIHUE.nam'
if platform.system() == 'Windows': exe_name='mfnwt.exe'
if platform.system() == 'Darwin': exe_name='mfnwt'


ml = flopy.modflow.Modflow.load(modelname, version='mfnwt', exe_name=exe_name, model_ws=MODFLOW_WORKSPACE)


coupled_model = CoupledModel(ml, 'llanquihue/SWMM/GIS/SWMM_S.shp')

gdf_final = coupled_model.joined_data

# gdf_final = georeference_models(ml, 'SWMM_inputs/shapes/SWMM_S.shp')
print(gdf_final.head())

gdf_final.plot(column="S", legend=True)
gdf_final.plot(column="SU", legend=True)
gdf_final.plot(column="subcatchment", legend=True)
plt.show()