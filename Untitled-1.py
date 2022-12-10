


import math
import os
import platform
import re
import sys
import time
from datetime import datetime

import fiona
import flopy
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flopy.discretization.structuredgrid import StructuredGrid
from flopy.utils import Raster
from geopandas import GeoDataFrame
from geopandas.tools import sjoin
from matplotlib.colors import LightSource
from numpy import *
from pandas import ExcelWriter
from pyproj import Proj, transform
from pyswmm import Links, Nodes, RainGages, Simulation, Subcatchments, SystemStats
from shapely.geometry import Point, Polygon

os.environ['KMP_DUPLICATE_LIB_OK']='True'

modelname='LLANQUIHUE'
if platform.system() == 'Windows': exe_name='mfnwt.exe'
if platform.system() == 'Darwin': exe_name='mfnwt'

workspace = 'LLANQUIHUE_MF_workspace'
if not os.path.isdir(workspace):
    os.mkdir(workspace)
    
ml = flopy.modflow.Modflow(modelname, version='mfnwt', exe_name=exe_name, model_ws=workspace)


nlay = 1     #number of layers
nrows = 100  #number of rows
ncols = 100  #number of columns
W = 5000     #model size (X)
H = 5000     #model size (Y)
d = 50      #aquifer thickness
top = 0
bottom = -50

x_resolution = W/ncols    
y_resolution = H/nrows

delr = y_resolution
delc = x_resolution


nper = 1
perlen = 1
nstp = 1
steady = True

dis = flopy.modflow.ModflowDis(
    ml, nlay=nlay, nrow=nrows, ncol=ncols,
    delr=delr, delc=delc, top=top, botm=bottom,
    nper=nper, perlen=perlen, nstp=nstp, steady=steady)


ml.modelgrid.set_coord_info(xoff=663776, yoff=5428040, angrot=0, epsg=32718)

ml.update_modelgrid()


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

mapview = flopy.plot.PlotMapView(model=ml)

linecollection = mapview.plot_grid()

t = ax.set_title("Model Grid")


ml.modelgrid.write_shapefile('GIS/SHAPES/MODFLOW_grid_flopy.shp')


MODFLOW_gdf=gpd.read_file('GIS/SHAPES/MODFLOW_grid_flopy.shp')

MODFLOW_gdf.crs


fig = plt.figure(figsize=(10, 6))

ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)


DEMraster_path='../GIS/DEM/clipped_dem.tif'

L1_top = np.ones((nrows,ncols)) 

for index, row in MODFLOW_gdf.iterrows():
    L1_top[row['row']-1,row['column']-1] = row['elevation']


dis = flopy.modflow.ModflowDis(
    ml, nlay=nlay, nrow=nrows, ncol=ncols,
    delr=delr, delc=delc, top=L1_top, botm=bottom,
    nper=nper, perlen=perlen, nstp=nstp, steady=steady)


dis.plot()

lakeBC_gdf=gpd.read_file('GIS/SHAPES/lake.shp')

riverBC_gdf=gpd.read_file('GIS/SHAPES/rio_maullin.shp')

riverBC_gdf.plot(), lakeBC_gdf.plot()

inactive_BC_gdf=gpd.read_file('GIS/SHAPES/inactive_cells.shp')

inactive_BC_gdf.plot()

MODFLOW_gdf['ibound']=1

poly_MODFLOWgrid_gdf = gpd.GeoDataFrame.from_file('GIS/SHAPES/MODFLOW_grid_flopy.shp')

poly_inactiveBC_gdf=gpd.GeoDataFrame.from_file('GIS/SHAPES/inactive_cells.shp')
poly_lakeBC_gdf=gpd.GeoDataFrame.from_file('GIS/SHAPES/lake.shp')
poly_riverBC_gdf=gpd.GeoDataFrame.from_file('GIS/SHAPES/rio_maullin.shp')

contained_inactivedBC_gdf=gpd.sjoin(poly_MODFLOWgrid_gdf, poly_inactiveBC_gdf, how="inner", op='intersects')
contained_riverBC_gdf=gpd.sjoin(poly_MODFLOWgrid_gdf, poly_riverBC_gdf, how="inner", op='intersects')
contained_lakeBC_gdf=gpd.sjoin(poly_MODFLOWgrid_gdf,poly_lakeBC_gdf, how="inner", op='intersects')

MODFLOW_gdf['inactive']=MODFLOW_gdf.index.isin(contained_inactivedBC_gdf.index.values.tolist())
MODFLOW_gdf['river']=MODFLOW_gdf.index.isin(contained_riverBC_gdf.index.values.tolist())
MODFLOW_gdf['lake']=MODFLOW_gdf.index.isin(contained_lakeBC_gdf.index.values.tolist())


MODFLOW_gdf.loc[MODFLOW_gdf.inactive == True, 'ibound'] = 0
MODFLOW_gdf.loc[MODFLOW_gdf.river == True, 'ibound'] = -1
MODFLOW_gdf.loc[MODFLOW_gdf.lake == True, 'ibound'] = -1

for i in range(len(MODFLOW_gdf["ibound"])):
    if MODFLOW_gdf["ibound"][i]==-1:
        MODFLOW_gdf["inactive"][i]=False

for i in range(len(MODFLOW_gdf["lake"])):
    if MODFLOW_gdf["lake"][i]:
        MODFLOW_gdf["river"][i]=False

plt.axis("off")


MODFLOW_gdf[['row', 'column']] = MODFLOW_gdf[['row', 'column']].astype(int)

ibound = np.ones((nlay,nrows,ncols)) 

for index, row in MODFLOW_gdf.iterrows():
    ibound[:,row['row']-1,row['column']-1] = row['ibound']

bas = flopy.modflow.ModflowBas(ml, ibound=ibound, strt=L1_top)

bas.plot()


upw = flopy.modflow.ModflowUpw(ml, hk=22.0, sy=0.15, laytyp=1, iphdry=0, ipakcb=53)


SWMM_S_gdf=gpd.read_file('GIS/SHAPES/SWMM_S.shp')

SWMM_SE_gdf=gpd.read_file('GIS/SHAPES/SWMM_SE.shp')
SWMM_SE_gdf["%imperv"]=0

imperv_S=pd.read_excel("S_imperv.xlsx")

SWMM_S_gdf=SWMM_S_gdf.merge(imperv_S, left_on='SWMM_NAME', right_on='NAME')

SWMM_S_gdf=SWMM_S_gdf.drop(columns=["NAME"])
SWMM_S_gdf=SWMM_S_gdf.rename(columns={"ID": "id"})

SWMM_S_gdf=SWMM_S_gdf.append(SWMM_SE_gdf)


MODFLOW_centroid=gpd.GeoDataFrame()
MODFLOW_centroid["geometry"]=MODFLOW_gdf["geometry"].centroid
MODFLOW_centroid["node"]=MODFLOW_gdf["node"]
MODFLOW_centroid.crs = {'init': 'epsg:32718'}

map_S_to_MFgrid_gdf=gpd.sjoin(MODFLOW_centroid, SWMM_S_gdf, how="left", op='intersects')
map_S_to_MFgrid_gdf.plot(column="%imperv", legend=True)


MODFLOW_gdf_SWMM=pd.DataFrame()
MODFLOW_gdf_SWMM=MODFLOW_gdf.join(map_S_to_MFgrid_gdf,rsuffix='_right')

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.drop(columns=["geometry_right","node_right", "index_right", "AREA"])
MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.rename(columns={"SWMM_NAME": "S"})

MODFLOW_gdf_SWMM.plot(column="id")


MODFLOW_gdf_SWMM["DRN"]=1

for i in range(len(MODFLOW_gdf_SWMM.groupby("id").count().index)):
    s=MODFLOW_gdf_SWMM.groupby("id").count().index[i]
    imperv=MODFLOW_gdf_SWMM[MODFLOW_gdf_SWMM["id"]==s]["%imperv"].mean()
    count_cells=MODFLOW_gdf_SWMM.groupby("id").count()["node"][s]
    imperv_cells=round(count_cells*imperv, 0)
    contador=0
    for j in MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM.id == s].index:
        if contador>=imperv_cells:
            break
        MODFLOW_gdf_SWMM["DRN"][j]=0
        contador+=1

is_drn = MODFLOW_gdf_SWMM[MODFLOW_gdf_SWMM['DRN'] == 1]


ax = SWMM_S_gdf.plot(figsize=(8,8), edgecolor='black', color="white");
MODFLOW_gdf_SWMM.plot(ax=ax,column='DRN', legend=True, alpha=0.8)
plt.axis("off")


DRN_burn_depth=0.0

# Global parameters needed to calculate drain conductance (see reference MODELMUSE DRN package pane)
DRN_L = x_resolution 
DRN_W = y_resolution
DRN_M = 1
DRN_K = 0.05 #m/dia      

DRN_C=DRN_K*DRN_L*DRN_W/DRN_M


DRN_stress_period_data=[]

for index,row in is_drn.iterrows():
    new_row=[0,row["row"]-1,row["column"]-1,row["elevation"]-DRN_burn_depth,DRN_C]
    DRN_stress_period_data.append(new_row)




drn = flopy.modflow.ModflowDrn(ml, stress_period_data=DRN_stress_period_data, ipakcb=53)



rch = flopy.modflow.ModflowRch(ml, nrchop=3, rech=0.001)


spd = {(0,0): ['save head', 'save budget'], (1,0): ['save head', 'save budget']}
oc = flopy.modflow.ModflowOc(ml, stress_period_data=spd) 


nwt = flopy.modflow.ModflowNwt(ml)

ml.get_package_list()


ml.check()

ml.write_input()
ml.run_model()


fname = os.path.join('LLANQUIHUE_MF_workspace', 'LLANQUIHUE.hds')
headfile = flopy.utils.HeadFile(fname, model=ml)
heads = headfile.get_data()
heads[heads==1.e+30] = np.nan            # fix masked data 
heads[heads==-999.99] = np.nan 

np.nanmin(heads[0]), np.nanmax(heads[0])

(48.875, 71.0)

strt = heads[0]
strt1 = heads[0]

fig = plt.figure(figsize=(10, 10))

levels = np.arange(0, 150, 1)

ax = fig.add_subplot(1, 2, 1, aspect='equal')
ax.set_title('Groundwater Heads-Initial Condition')
mapview = flopy.plot.PlotMapView(model=ml)
quadmesh = mapview.plot_ibound()
quadmesh = mapview.plot_array(heads, masked_values=[999.], alpha=0.5)
plt.axis("off")


boundry_conditions_df=pd.read_excel("elev_maullin_llanq.xlsx", index_col=0)
boundry_conditions_df.plot(figsize=(10,8))

dire="SWMM_Llanquihue"
file="/Llanquihue_base.inp"
SWMM_path=dire + "/" + file

with Simulation(SWMM_path) as sim:
    if not sim.flow_units=="CMS":
        print("SWMM flow units must be CMS")


SWMM_names_df=pd.DataFrame()

SWMM_names_df=pd.read_excel("SWMM_elem_names.xlsx")

#SWMM Subcatchments names:
S_names_list=SWMM_names_df["S_names"].dropna(how=all).astype(str).values.tolist()
#SWMM Wetlands Storge Unites names:
WSU_names_list=SWMM_names_df["WSU_names"].dropna(how=all).astype(str).values.tolist()
#Wetlands subcatchments
SW_names_list=SWMM_names_df["SW_names"].dropna(how=all).astype(str).values.tolist()
#SWMM Nodes names:
J_names_list=SWMM_names_df["J_names"].dropna(how=all).astype(str).values.tolist()
#SWMM conduits names:
C_names_list=SWMM_names_df["C_names"].dropna(how=all).astype(str).values.tolist()
#SWMM weirs names:
V_names_list=SWMM_names_df["V_names"].dropna(how=all).astype(str).values.tolist()
#SWMM SU names:
FSU_names_list=SWMM_names_df["FSU_names"].dropna(how=all).astype(str).values.tolist()


with Simulation(SWMM_path) as sim:
    S_list=[]
    WSU_list=[]
    FSU_list=[]
    J_list=[]
    C_list=[]
    SW_list=[]
    V_list=[]
    S_areas=[]
    for s in S_names_list:
        S_list.append(Subcatchments(sim)[s])
        S_areas.append(Subcatchments(sim)[s].area) #[ha]  
    for i in range(len(WSU_names_list)):
        WSU_list.append(Nodes(sim)[WSU_names_list[i]])
    for i in range(len(FSU_names_list)):
        FSU_list.append(Nodes(sim)[FSU_names_list[i]])
    for j in J_names_list:
        J_list.append(Nodes(sim)[j])
    for s in SW_names_list:
        SW_list.append(Subcatchments(sim)[s])
    for c in C_names_list:
        C_list.append(Links(sim)[c])
    for v in V_names_list:
        V_list.append(Links(sim)[v])


S_SWMM_TS=[]
WSU_SWMM_TS=[]
FSU_SWMM_TS=[]
SW_SWMM_TS=[]
J_SWMM_TS=[]
C_SWMM_TS=[]
V_SWMM_TS=[]
for i in S_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Precipitation (mm/h)"]=""
    a["Evaporation (mm/d)"]=""
    a["Infiltration (mm/h)"]="" 
    a["Runoff (m3/s)"]=""
    a["Runon (m3/s)"]=""
    a["Cumulative Infiltration (m3)"]=""
    a["Cumulative Evaporation (m3)"]=""
    S_SWMM_TS.append(a)

for i in WSU_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Head (m)"]=""
    a["Flooding (m3/s)"]=""
    a["Lateral inflow (m3/s)"]=""
    a["Total inflow (m3/s)"]=""
    a["Total outflow (m3/s)"]=""
    a["Volume (m3)"]=""
    a["Losses (m3/s)"]=""
    a["Cumulative Exfiltration Loss (m3)"]=""
    a["Cumulative Evaporation Loss (m3)"]=""
    WSU_SWMM_TS.append(a)
for i in FSU_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Head (m)"]=""
    a["Flooding (m3/s)"]=""
    a["Lateral inflow (m3/s)"]=""
    a["Total inflow (m3/s)"]=""
    a["Total outflow (m3/s)"]=""
    a["Volume (m3)"]=""
    a["Losses (m3/s)"]=""
    a["Cumulative Exfiltration Loss (m3)"]=""
    a["Cumulative Evaporation Loss (m3)"]=""
    FSU_SWMM_TS.append(a)
for i in J_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Head (m)"]=""
    a["Flooding (m3/s)"]=""
    a["Lateral inflow (m3/s)"]=""
    a["Total inflow (m3/s)"]=""
    a["Total outflow (m3/s)"]=""
    J_SWMM_TS.append(a)
for i in C_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Flow (m3/s)"]=""
    C_SWMM_TS.append(a)
for i in V_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Flow (m3/s)"]=""
    V_SWMM_TS.append(a)


MODFLOW_gdf_SWMM.plot(column="id")


extern_active_gdf=gpd.read_file('GIS/SHAPES/extern_active.shp')
extern_active_gdf.plot()


map_extern_to_MFgrid_gdf=gpd.sjoin(MODFLOW_centroid, extern_active_gdf, how="left", op='intersects')

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.join(map_extern_to_MFgrid_gdf,rsuffix='_right')

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.drop(columns=["geometry_right","node_right","index_right"])             

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.rename(columns={"id_right": "active_extern"})
MODFLOW_gdf_SWMM.plot(column="active_extern")



MODFLOW_gdf_SWMM["Nan_cell"]=0
MODFLOW_gdf_SWMM["Double_cell"]=0
for i in range(len(MODFLOW_gdf_SWMM["DRN"])):
    if MODFLOW_gdf_SWMM["ibound"][i]==1 and type(MODFLOW_gdf_SWMM["S"][i])!=str and MODFLOW_gdf_SWMM["active_extern"][i]!=1:
        MODFLOW_gdf_SWMM["Nan_cell"][i]=1
    if MODFLOW_gdf_SWMM["ibound"][i]==-1 and type(MODFLOW_gdf_SWMM["S"][i])==str:
        MODFLOW_gdf_SWMM["Double_cell"][i]=1
    if MODFLOW_gdf_SWMM["active_extern"][i]==1 and type(MODFLOW_gdf_SWMM["S"][i])==str:
        MODFLOW_gdf_SWMM["Double_cell"][i]=1
    if MODFLOW_gdf_SWMM["river"][i] and MODFLOW_gdf_SWMM["lake"][i]:
        MODFLOW_gdf_SWMM["Double_cell"][i]=1
    if (MODFLOW_gdf_SWMM["river"][i] or MODFLOW_gdf_SWMM["lake"][i]) and MODFLOW_gdf_SWMM["inactive"][i]:
        MODFLOW_gdf_SWMM["Double_cell"][i]=1

for i in range(len(MODFLOW_gdf_SWMM["DRN"])):
    if MODFLOW_gdf_SWMM["ibound"][i]==1 and type(MODFLOW_gdf_SWMM["S"][i])!=str and MODFLOW_gdf_SWMM["active_extern"][i]!=1:
        MODFLOW_gdf_SWMM["S"][i]=MODFLOW_gdf_SWMM["S"][i-1]
        MODFLOW_gdf_SWMM["id"][i]=MODFLOW_gdf_SWMM["id"][i-1]    
    if MODFLOW_gdf_SWMM["ibound"][i]==-1 and type(MODFLOW_gdf_SWMM["S"][i])==str:
        MODFLOW_gdf_SWMM["S"][i]=np.nan
        MODFLOW_gdf_SWMM["id"][i]=np.nan

S_areas_cells=MODFLOW_gdf_SWMM.groupby("S").node.count()*x_resolution*y_resolution

S_areas_modflow=[]
for s in S_names_list:
    if s in S_areas_cells.index:
        S_areas_modflow.append(S_areas_cells[s])

for i in range(len(S_areas)):
    print(S_names_list[i], S_areas[i]*10000/S_areas_modflow[i])

MODFLOW_gdf_SWMM["WSU"]=""
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM["S"] == "S_W1", 'WSU'] = "W1"
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM["S"] == "S_W2", 'WSU'] = "W2"
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM["S"] == "S_W3", 'WSU'] = "W3"
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM["S"] == "S_W4", 'WSU'] = "W4"
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM["S"] == "S_W5", 'WSU'] = "W5"

WSU_areas_cells=MODFLOW_gdf_SWMM.groupby("WSU").node.count()*x_resolution*y_resolution
WSU_areas_modflow=[]
for wsu in WSU_names_list:
    if wsu in WSU_areas_cells.index:
        WSU_areas_modflow.append(WSU_areas_cells[wsu])


FSU_gdf=gpd.read_file('GIS/SHAPES/flood_zone.shp')
FSU_gdf.plot()


map_FSU_to_MFgrid_gdf=gpd.sjoin(MODFLOW_centroid,FSU_gdf, how="left", op='intersects')
map_FSU_to_MFgrid_gdf.plot(column="id")

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.join(map_FSU_to_MFgrid_gdf,rsuffix='_right')

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.drop(columns=["geometry_right","node_right","index_right","id_right"])      
MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.rename(columns={"SWMM_NAME": "FSU"})

FSU_areas_cells=MODFLOW_gdf_SWMM.groupby("FSU").node.count()*x_resolution*y_resolution
FSU_areas_modflow=[]
for fsu in FSU_names_list:
    if fsu in FSU_areas_cells.index:
        FSU_areas_modflow.append(FSU_areas_cells[fsu])

SWMM_TS1_gdf=gpd.read_file('GIS/SHAPES/Teodosio_Sarao_zone1.shp')
SWMM_TS2_gdf=gpd.read_file('GIS/SHAPES/Teodosio_Sarao_zone2.shp')
SWMM_TS1_gdf.plot()
SWMM_TS2_gdf.plot()


map_TS1_to_MFcentroid_gdf=gpd.sjoin(MODFLOW_centroid, SWMM_TS1_gdf, how="inner", op='intersects')
map_TS1_to_MFcentroid_gdf.plot()

map_TS2_to_MFcentroid_gdf=gpd.sjoin(MODFLOW_centroid, SWMM_TS2_gdf, how="inner", op='intersects')
map_TS2_to_MFcentroid_gdf.plot()


SWMM_pointsTS_gdf = gpd.read_file('GIS/SHAPES/points_teodosio_sarao.shp')
SWMM_pointsTS_gdf.plot()


gpd1 = map_TS1_to_MFcentroid_gdf
gpd2 = map_TS2_to_MFcentroid_gdf
gpd3 = SWMM_pointsTS_gdf

def ckdnearest(gdA, gdB):
    nA = np.array(list(zip(gdA.centroid.x, gdA.centroid.y)) )
    nB = np.array(list(zip(gdB.geometry.x, gdB.geometry.y)) )
    btree = cKDTree(nB)
    dist, idx = btree.query(nA, k=1)
    gdf = pd.concat(
        [gdA.reset_index(drop=True), gdB.loc[idx, gdB.columns != 'geometry'].reset_index(drop=True),
         pd.Series(dist, name='dist')], axis=1)
    return gdf

map_MFdraincells_to_TS1routingnodes_gdf = ckdnearest(gpd1, gpd3)
map_MFdraincells_to_TS2routingnodes_gdf = ckdnearest(gpd2, gpd3)

map_MFdraincells_to_TS1routingnodes_gdf.plot(column='fid')
map_MFdraincells_to_TS2routingnodes_gdf.plot(column='fid')

map_MFdraincells_to_TS_gdf=pd.concat([map_MFdraincells_to_TS1routingnodes_gdf, map_MFdraincells_to_TS2routingnodes_gdf])
map_MFdraincells_to_TS_gdf=map_MFdraincells_to_TS_gdf.drop(["fid", "geometry","index_right","ID","dist", "AREA", "SWMM_NAME"], axis=1)
map_MFdraincells_to_TS_gdf=map_MFdraincells_to_TS_gdf.rename(columns={"SWMM_NAM_1": "drn_to_TS"})
map_MFdraincells_to_TS_gdf

drnto_df = pd.read_excel("drn_to.xlsx")

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.merge(drnto_df, on="S", how="left" )

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.merge(map_MFdraincells_to_TS_gdf, on="node", how="left" )

contador=0
for i in range(len(MODFLOW_gdf_SWMM.drn_to)):
    if type(MODFLOW_gdf_SWMM.drn_to[i])!=str:
        if type(MODFLOW_gdf_SWMM.drn_to_TS[i])==str:
            MODFLOW_gdf_SWMM.drn_to[i]=MODFLOW_gdf_SWMM.drn_to_TS[i]
            contador+=1
print(contador)

MODFLOW_gdf_SWMM=MODFLOW_gdf_SWMM.drop(columns=["drn_to_TS"])

contador=0
for i in range(len(MODFLOW_gdf_SWMM.drn_to)):
    if MODFLOW_gdf_SWMM["ibound"][i]==1:
        if not (MODFLOW_gdf_SWMM["active_extern"][i]==1):
            if MODFLOW_gdf_SWMM["DRN"][i]==1 and type(MODFLOW_gdf_SWMM["drn_to"][i])!=str:
                contador+=1
                print(i)

MODFLOW_gdf_SWMM["zon"]=1
aliases = {1: 'Active_Extern' }

contador=2
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM.lake == True, 'zon'] = contador
aliases.update({contador:"Lake"})

contador=contador+1
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM.river == True, 'zon'] = contador
aliases.update({contador:"River"})

contador=contador+1
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM.inactive == True, 'zon'] = contador
aliases.update({contador:"Inactive"})

contador=contador+1
for i in range(len(S_names_list)):
    MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM.S == S_names_list[i], 'zon'] = contador
    aliases.update({contador: S_names_list[i]})
    contador=contador+1

MODFLOW_gdf_SWMM.plot(column="zon", legend=True)

df_aliases = pd.DataFrame.from_dict(aliases, orient='index')
df_aliases.to_csv('aliases.csv', index=False)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

ax = MODFLOW_gdf_SWMM.plot(ax=ax, column="zon",  cmap='Spectral', alpha=0.9)
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
#plt.colorbar(ax.images[0], shrink=0.7)
pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
pmv.plot_grid(ax=ax, alpha=0.2)
plt.axis("off")

zon = np.zeros((nlay,nrows,ncols)) 

for index, row in MODFLOW_gdf_SWMM.iterrows():
    zon[:,row['row']-1,row['column']-1] = row['zon']

zon=zon[0].astype(np.int64)

MODFLOW_gdf_SWMM["drn_to2"]=np.nan
MODFLOW_gdf_SWMM.loc[MODFLOW_gdf_SWMM.DRN == 1, 'drn_to2'] = MODFLOW_gdf_SWMM.drn_to

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')


#ax = SWMM_S_gdf.plot(ax=ax,figsize=(10,12), edgecolor='black', color="white",alpha=0.8);
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
ax = MODFLOW_gdf_SWMM.plot(ax=ax, column="drn_to2",  cmap='Spectral')
#plt.colorbar(ax.images[0], shrink=0.7)
pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
pmv.plot_grid(ax=ax, alpha=0.2);
plt.axis("off")

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

ax = MODFLOW_gdf.plot(ax=ax, column="elevation",  cmap='Spectral', legend=False, vmin=45, vmax=145)
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
pmv.plot_grid(ax=ax, alpha=0.2);
plt.axis("off")
norm = colors.Normalize(vmin=45, vmax=145)
cbar = plt.cm.ScalarMappable(norm=norm, cmap='Spectral')
ax_cbar = fig.colorbar(cbar, ax=ax, label="Elevation (m.a.s.l)")

MODFLOW_gdf_SWMM["strt"]=""
strt_df=np.reshape(strt1, len(MODFLOW_gdf_SWMM))
MODFLOW_gdf_SWMM["strt"]=strt_df


for i in range(len(MODFLOW_gdf_SWMM.strt)):
    if MODFLOW_gdf_SWMM["river"][i]:
        MODFLOW_gdf_SWMM["strt"][i]=boundry_conditions_df.Elev_maullin["2017-03-01 00:00"]
    if MODFLOW_gdf_SWMM["lake"][i]:
        MODFLOW_gdf_SWMM["strt"][i]=boundry_conditions_df.Elev_llanq["2017-03-01 00:00"]


MODFLOW_gdf_SWMM.plot(column="strt", legend=True)
plt.axis("off")

perlen,nstp,nper,steady = 1, 1, 1, False
dis = flopy.modflow.ModflowDis(ml, nlay=nlay, nrow=nrows, ncol=ncols,delr=delr, delc=delc, top=L1_top, botm=bottom, nper=nper, perlen=perlen, nstp=nstp, steady=steady);


ZB_TS=[]

TS_gdf=[]

with Simulation(SWMM_path) as sim:
    
    system_routing = SystemStats(sim)
    #Lists of cummulative infiltration to calculate delta infiltratation for every time step: S, WSU and FSU 
    inf_S_list_1=np.zeros(len(S_list))
    inf_S_list_2=np.zeros(len(S_list))
    inf_WSU_list_1=np.zeros(len(WSU_list))
    inf_WSU_list_2=np.zeros(len(WSU_list))
    inf_FSU_list_1=np.zeros(len(FSU_list))
    inf_FSU_list_2=np.zeros(len(FSU_list))
    
    #Lists for the DRN incorporation
    rate_J_list=np.zeros(len(J_list))
    rate_WSU_list=np.zeros(len(WSU_list))
    rate_FSU_list=np.zeros(len(FSU_list))
    
    #time counter for daily infiltration agregation and hourly reports
    step_counter=0
    day_counter=0
    hourly_counter=0
    for step in sim:
        step_counter=step_counter+1

        if step_counter==360: #CHANGE ACCORDING TO DT
            step_counter=0
            hourly_counter+=1
            
            # TIME SERIES RESULTS: 1HR AS REPORT TIME STEP
            
            #SUBCATCHMENTS TIME SERIES RESULTS
            for i in range(len(S_list)):
                new_row = {'Time': sim.current_time, "Precipitation (mm/h)":S_list[i].rainfall, 
                           "Evaporation (mm/d)":S_list[i].evaporation_loss,"Infiltration (mm/h)":S_list[i].infiltration_loss, 
                           "Runoff (m3/s)":S_list[i].runoff,"Runon (m3/s)":S_list[i].runon,
                           "Cumulative Infiltration (m3)": S_list[i].statistics["infiltration"], "Cumulative Evaporation (m3)": S_list[i].statistics["evaporation"]}
                S_SWMM_TS[i] =S_SWMM_TS[i].append(new_row, ignore_index=True)

            #WETLAND STORAGE UNITS TIME SERIES RESULTS
            for i in range(len(WSU_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":WSU_list[i].depth, 
                           "Head (m)":WSU_list[i].head, "Flooding (m3/s)":WSU_list[i].flooding, 
                           "Lateral inflow (m3/s)":WSU_list[i].lateral_inflow,"Total inflow (m3/s)":WSU_list[i].total_inflow,
                           "Total outflow (m3/s)":WSU_list[i].total_outflow, "Volume (m3)":WSU_list[i].volume, "Losses (m3/s)":WSU_list[i].losses,
                           "Cumulative Exfiltration Loss (m3)": WSU_list[i].storage_statistics["exfil_loss"], "Cumulative Evaporation Loss (m3)": WSU_list[i].storage_statistics["evap_loss"]}
                WSU_SWMM_TS[i] =WSU_SWMM_TS[i].append(new_row, ignore_index=True)
            #FLOOD STORAGE UNITS TIME SERIES RESULTS
            for i in range(len(FSU_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":FSU_list[i].depth, 
                           "Head (m)":FSU_list[i].head, "Flooding (m3/s)":FSU_list[i].flooding, 
                           "Lateral inflow (m3/s)":FSU_list[i].lateral_inflow,"Total inflow (m3/s)":FSU_list[i].total_inflow,
                           "Total outflow (m3/s)":FSU_list[i].total_outflow, "Volume (m3)":FSU_list[i].volume, "Losses (m3/s)":FSU_list[i].losses,
                           "Cumulative Exfiltration Loss (m3)": FSU_list[i].storage_statistics["exfil_loss"], "Cumulative Evaporation Loss (m3)": FSU_list[i].storage_statistics["evap_loss"]}
                FSU_SWMM_TS[i] =FSU_SWMM_TS[i].append(new_row, ignore_index=True)
                
            #NODES TIME SERIES RESULTS
            for i in range(len(J_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":J_list[i].depth, 
                          "Head (m)":J_list[i].head, "Flooding (m3/s)":J_list[i].flooding, 
                          "Lateral inflow (m3/s)":J_list[i].lateral_inflow,"Total inflow (m3/s)":J_list[i].total_inflow,
                          "Total outflow (m3/s)":J_list[i].total_outflow}
                J_SWMM_TS[i] =J_SWMM_TS[i].append(new_row, ignore_index=True) 
                
                
            #CONDUITS TIME SERIES RESULTS
            for i in range(len(C_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":C_list[i].depth, 
                           "Flow (m3/s)":C_list[i].flow}
                C_SWMM_TS[i] =C_SWMM_TS[i].append(new_row, ignore_index=True)
            
            # WEIRS TIME SERIES RESULTS
            for i in range(len(V_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":V_list[i].depth, 
                          "Flow (m3/s)":V_list[i].flow}
                V_SWMM_TS[i] =V_SWMM_TS[i].append(new_row, ignore_index=True)          
        
        #DAILY INFILTRATION ON S, WSU AND FSU:
        
        if hourly_counter==24:
            day_counter=day_counter+1
            hourly_counter=0
            
            print(sim.current_time)

            for i in range(len(S_list)):
                #Delta infiltration
                
                inf_S_list_2[i]=(S_list[i].statistics["infiltration"]-inf_S_list_1[i])
                inf_S_list_1[i]=S_list[i].statistics["infiltration"]
                
            for i in range(len(WSU_list)):
                #Delta infiltration
                inf_WSU_list_2[i]=(WSU_list[i].storage_statistics["exfil_loss"]-inf_WSU_list_1[i])
                inf_WSU_list_1[i]=WSU_list[i].storage_statistics["exfil_loss"]

            for i in range(len(FSU_list)):
                #Delta infiltration
                inf_FSU_list_2[i]=(FSU_list[i].storage_statistics["exfil_loss"]-inf_FSU_list_1[i])
                inf_FSU_list_1[i]=FSU_list[i].storage_statistics["exfil_loss"]

                
            RCH_S=inf_S_list_2
            RCH_WSU=inf_WSU_list_2
            RCH_FSU=inf_FSU_list_2

            #CHANGE OF UNITS m3/day->m/day:
    
            RCH_S_M=RCH_S/(np.array(S_areas_modflow))
            RCH_WSU_M=RCH_WSU/(np.array(WSU_areas_modflow))
            RCH_FSU_M=RCH_FSU/(np.array(FSU_areas_modflow))            
        
            RCH_S_df=pd.DataFrame({"S":S_names_list, "RCH_S":RCH_S_M})
            RCH_WSU_df=pd.DataFrame({"WSU":WSU_names_list, "RCH_WSU":RCH_WSU_M})
            RCH_FSU_df=pd.DataFrame({"FSU":FSU_names_list, "RCH_FSU":RCH_FSU_M})
            
            #Gereferenced RCH: Add to the MODFLOW_gdf_test new columns RCH_S, RCH_WSU, RCH_FSU
            
            MODFLOW_gdf_loop=pd.DataFrame()
            
            MODFLOW_gdf_loop=pd.merge(MODFLOW_gdf_SWMM, RCH_S_df, on="S", how="left")
            MODFLOW_gdf_loop=pd.merge(MODFLOW_gdf_loop, RCH_WSU_df, on="WSU", how="left")
            MODFLOW_gdf_loop=pd.merge(MODFLOW_gdf_loop, RCH_FSU_df, on="FSU", how="left")
            
            # Sum georeferences RCHs
            
            MODFLOW_gdf_loop["RCH"]= MODFLOW_gdf_loop.RCH_S.fillna(0) + MODFLOW_gdf_loop.RCH_WSU.fillna(0) +  MODFLOW_gdf_loop.RCH_FSU.fillna(0)
            
            #Create MODFLOW inputs: RCH package
            
            rch_array = np.zeros((nrows,ncols))
            recharge_cells = MODFLOW_gdf_SWMM.index.values
            for cell in recharge_cells:
                row = MODFLOW_gdf_loop.row[cell]
                col = MODFLOW_gdf_loop.column[cell]
                flux = MODFLOW_gdf_loop.RCH[cell]
                rch_array[row - 1][col - 1] = flux
            rch_array[np.isnan(rch_array)] = 0.
        
            rch = flopy.modflow.ModflowRch(ml, nrchop=3,rech=rch_array, ipakcb=53)
            
            #Update boundry conditions for lake and river
            
            MODFLOW_gdf_loop["strt"]=""
            strt_df=np.reshape(strt, len(MODFLOW_gdf_loop))
            MODFLOW_gdf_loop["strt"]=strt_df
            for i in range(len(MODFLOW_gdf_loop.strt)):
                if MODFLOW_gdf_loop["river"][i]:
                    MODFLOW_gdf_loop["strt"][i]=boundry_conditions_df.Elev_maullin[sim.current_time]
                if MODFLOW_gdf_loop["lake"][i]:
                    MODFLOW_gdf_loop["strt"][i]=boundry_conditions_df.Elev_llanq[sim.current_time]
                    
            strt_array = np.zeros((nrows,ncols))
            strt_cells = MODFLOW_gdf_SWMM.index.values
            for cell in strt_cells:
                row = MODFLOW_gdf_loop.row[cell]
                col = MODFLOW_gdf_loop.column[cell]
                flux = MODFLOW_gdf_loop.strt[cell]
                strt_array[row - 1][col - 1] = flux
            
            #bas = flopy.modflow.ModflowBas(ml, ibound=ibound, strt=strt)
            bas = flopy.modflow.ModflowBas(ml, ibound=ibound, strt=strt_array) #use the head table of the last time step and bc
           
            #Run MODFLOW
            
            ml.write_input()
            ml.run_model(silent=True)
            
            #Read MODFLOW outputs
            fname = os.path.join('LLANQUIHUE_MF_workspace', 'LLANQUIHUE.hds')
            headfile = flopy.utils.HeadFile(fname, model=ml)
            heads = headfile.get_data()
            heads[heads==1.e+30] = np.nan            # fix masked data 
            heads[heads==-999.99] = np.nan

            #Strt next loop
            
            strt = heads[0]
            
            top = ml.dis.top.array
            DTWT = (top-DRN_burn_depth)-heads[0]  
            
            #DRN calculation
            
            delta_H=np.reshape(DTWT, len(MODFLOW_gdf_loop))
            altura=np.reshape(heads[0], len(MODFLOW_gdf_loop))
            for i in range(len(delta_H)):
                if delta_H[i]<0:
                    delta_H[i]=-delta_H[i]
                else:
                     delta_H[i]=0.
            MODFLOW_gdf_loop["Altura"]=altura   
            MODFLOW_gdf_loop["delta_H"]=delta_H
            MODFLOW_gdf_loop["DRN_rate"]=0.
            
            for i in range(len(MODFLOW_gdf_loop["DRN"])):
                if MODFLOW_gdf_loop["DRN"][i]==1:
                    MODFLOW_gdf_loop["DRN_rate"][i]=(delta_H[i])*DRN_C
                if MODFLOW_gdf_loop["ibound"][i]==-1:
                    MODFLOW_gdf_loop["DRN_rate"][i]=0
                
            
            #INFLOW RATES IN WSU, FSU AND JUNCTIONS
            
            #Inflow rates in WSU:
            
            inflow_WSU=MODFLOW_gdf_loop.groupby("drn_to").sum()["DRN_rate"]
            inflow_WSU_list=[]
            for w in WSU_names_list:
                if w in inflow_WSU.index:
                    inflow_WSU_list.append(inflow_WSU[w])
                else:
                    inflow_WSU_list.append(0.)
                    
            #Inflow rates in FSU:
            
            inflow_FSU=MODFLOW_gdf_loop.groupby("drn_to").sum()["DRN_rate"]
            inflow_FSU_list=[]
            for f in FSU_names_list:
                if f in inflow_FSU.index:
                    inflow_FSU_list.append(inflow_FSU[f])
                else:
                    inflow_WSU_list.append(0.)
                    
            #Inflow rate in Junctions
            
            inflow_j=MODFLOW_gdf_loop.groupby("drn_to").sum()["DRN_rate"]
            inflow_j_list=[]
            for j in J_names_list:
                if j in inflow_j.index:
                    inflow_j_list.append(inflow_j[j])
                else:
                    inflow_j_list.append(0.)
                
            #Generated inflow SWMM WSU:
            
            for i in range(len(WSU_list)):
                rate=inflow_WSU_list[i]/86400. #m3/s
                WSU_list[i].generated_inflow(rate)
                
            #Generated inflow SWMM FSU:
            
            for i in range(len(FSU_list)):
                rate=inflow_FSU_list[i]/86400. #m3/s
                FSU_list[i].generated_inflow(rate)
                
            #Generated inflow SWMM Junctions:
            
            for i in range(len(J_list)):
                rate=inflow_j_list[i]/86400 #m3/s
                J_list[i].generated_inflow(rate)
                
            #Save MODFLOW cells information
            TS_gdf.append(MODFLOW_gdf_loop)
            
                
            #Save zonebudget DataFrame
            fname = os.path.join('LLANQUIHUE_MF_workspace', 'LLANQUIHUE.cbc')
            cbb = flopy.utils.CellBudgetFile(fname)
            zb = flopy.utils.ZoneBudget(cbb, zon, aliases=aliases)
            zb_df=zb.get_dataframes()
            ZB_TS.append(zb_df)
            
                
    routing_stats=system_routing.routing_stats
    runoff_stats=system_routing.runoff_stats
            
print("Flow Routing Mass Balance Error:", sim.flow_routing_error)
print("Runoff Mass Balance Error:", sim.runoff_error)

Wall time: 4h 52min

Results

df_routing_couple = pd.DataFrame({'Dry Weather Inflow (m3)': [routing_stats['dry_weather_inflow']],
                   'Wet Weather Inflow (m3)': [routing_stats['wet_weather_inflow']], 
                   'Ground Water Inflow (m3)': [routing_stats['groundwater_inflow']],
                   #'RDII Inflow (m3)': [routing_stats['II_inflow']],
                   'External Inflow (m3)': [routing_stats['external_inflow']],
                   'External Outflow (m3)': [routing_stats['outflow']],
                   'Flooding Loss (m3)': [routing_stats['flooding']],
                   'Evaporation Loss (m3)': [routing_stats['evaporation_loss']],
                   'Seepage Loss (m3)': [routing_stats['seepage_loss']],
                   'Initial Stored Volume': [routing_stats['initial_storage']],
                   'Final Stored Volume': [routing_stats['final_storage']]})

df_runoff_couple = pd.DataFrame({'Total Precipitation (mm)': [runoff_stats['rainfall']],
                   'Evaporation Loss (m3)': [runoff_stats['evaporation']], 
                   'Infiltration Loss (m3)': [runoff_stats['infiltration']],
                   'Runoff (m3)': [runoff_stats['runoff']],
                   "Initial Storage (mm)":[runoff_stats['init_storage']],
                   "Final Storage (mm)":[runoff_stats['final_storage']]})

    Report for SWMM model

df_routing_couple

	Dry Weather Inflow (m3) 	Wet Weather Inflow (m3) 	Ground Water Inflow (m3) 	External Inflow (m3) 	External Outflow (m3) 	Flooding Loss (m3) 	Evaporation Loss (m3) 	Seepage Loss (m3) 	Initial Stored Volume 	Final Stored Volume
0 	0.0 	2.531414e+06 	0.0 	9.426430e+06 	6.635352e+06 	0.0 	70819.796037 	5.282097e+06 	2.668431e+06 	0.0

df_runoff_couple

	Total Precipitation (mm) 	Evaporation Loss (m3) 	Infiltration Loss (m3) 	Runoff (m3) 	Initial Storage (mm) 	Final Storage (mm)
0 	3228.453076 	1.303372e+06 	1.316862e+07 	2.531414e+06 	0.0 	0.0

    Working with ZB output:

ZB_TS[0]=ZB_TS[0].reset_index(level=[0,1])

ZB by zone follows aliases index

dias=np.linspace(0,int(len(ZB_TS)-1),int(len(ZB_TS)))
ZB_by_zone=[]
columnas=ZB_TS[0].name

dias=np.linspace(0,int(len(ZB_TS)-1),int(len(ZB_TS)))
ZB_by_zone=[]
columnas=ZB_TS[0].name
for z in aliases:
    zone=aliases[z]
    df_ZB=pd.DataFrame(index=dias)
    for i in columnas:
        df_ZB[i]=""
    for i in range(len(ZB_TS)):
        for j in range(len(columnas)):
            df_ZB[columnas[j]][i]=ZB_TS[i].iloc[j][zone]
    ZB_by_zone.append(df_ZB)

for i in range(len(ZB_by_zone)):
    ZB_by_zone[i]["PERCENT_DISCREPANCY"].plot(figsize=(10,3), label=aliases[i+1])
plt.ylabel("% Discrepancy")
plt.xlabel("Días")
#plt.legend()
plt.show()

Continuity Analysis
Coupling Continuity

    External Inflow in SWMM

E_in_s=df_routing_couple["External Inflow (m3)"][0]
E_in_s

9426429.608449133

    Drain water by MODFLOW using ZB:

DRN_S_m=0
contador=0
for i in range(len(aliases)):
    if aliases[i+1] in S_names_list:
        contador+=1
        DRN_S_m+=sum(list(ZB_by_zone[i]["TO_DRAINS"])[0:-1])

DRN_S_m

9426425.71767044

error_drn_inflow=((DRN_S_m)-E_in_s)/(DRN_S_m)*100
print("Continuity Coupling Error(DRN/INFLOW):",error_drn_inflow,"%")

Continuity Coupling Error(DRN/INFLOW): -4.127522784320767e-05 %

    Infiltration and Seepage Loss using SWMM output:

I_s=df_runoff_couple["Infiltration Loss (m3)"][0]
I_s

13168620.891738793

S_s=df_routing_couple["Seepage Loss (m3)"][0]
S_s

5282096.964601642

    Recharge in MODFLOW using ZB:

RCH_m=0
contador=0
for i in range(len(aliases)):
    if aliases[i+1] in S_names_list:
        contador+=1
        RCH_m+=sum(list(ZB_by_zone[i]["FROM_RECHARGE"]))

contador

47

error_inf_rch=(RCH_m-I_s-S_s)/RCH_m*100
print("Continuity Coupling Error(INF/RCH):",error_inf_rch,"%")

Continuity Coupling Error(INF/RCH): -1.5183199059592918e-05 %

SWMM Continuity

P_s=df_runoff_couple["Total Precipitation (mm)"][0]*sum(S_areas)*10
P_s

16998380.10915186

E_out_s=df_routing_couple["External Outflow (m3)"][0]
E_out_s

6635352.437175099

E_out_s=df_routing_couple["External Outflow (m3)"][0]
E_out_s

6635352.437175099

F_s=df_routing_couple["Flooding Loss (m3)"][0]
F_s

0.0

EVT_s=df_runoff_couple["Evaporation Loss (m3)"][0]+df_routing_couple["Evaporation Loss (m3)"][0]
EVT_s

1374191.4113278678

S_f_s= 44682+0*sum(S_areas)*10
S_i_s=75562
S_f_s

44682.0

SWMM_error=P_s+E_in_s-E_out_s-EVT_s-F_s-I_s-S_s-S_f_s+S_i_s
total_in_s=P_s+E_in_s
percent_SWMM_error=SWMM_error/total_in_s
print("SWMM Error:",percent_SWMM_error, "%" )

SWMM Error: -0.00017301873849875733 %

MODFLOW Continuity

CH_in_m=0
CH_out_m=0
S_to_m=0
S_from_m=0
DRN_m=0
for i in range(len(aliases)):
        CH_in_m+=sum(list(ZB_by_zone[i]["FROM_CONSTANT_HEAD"]))
        CH_out_m+=sum(list(ZB_by_zone[i]["TO_CONSTANT_HEAD"]))
        S_to_m+=sum(list(ZB_by_zone[i]["TO_STORAGE"]))
        S_from_m+=sum(list(ZB_by_zone[i]["FROM_STORAGE"]))
        DRN_m+=sum(list(ZB_by_zone[i]["TO_DRAINS"]))

MODFLOW_error=RCH_m+CH_in_m-DRN_m-CH_out_m-S_to_m+S_from_m
total_in_m=RCH_m+CH_in_m
percent_MODFLOW_error=MODFLOW_error/total_in_m*100
print("MODFLOW Error:",percent_MODFLOW_error, "%" )

MODFLOW Error: 0.22464729790435117 %

Global Continuity

DRN_ext_m=sum(list(ZB_by_zone[0]["TO_DRAINS"]))
DRN_ext_m

73508.67

dt_drn_m=0
for i in range(len(aliases)):
    if aliases[i+1] in S_names_list:
        dt_drn_m+=sum(list(ZB_by_zone[i]["TO_DRAINS"])[-1])
dt_drn_m        

6314.3847007751465

DRN_m-(dt_drn_m+DRN_ext_m+DRN_S_m)

0.2715492248535156

DRN_S_m-E_in_s

-3.8907786924391985

RCH_m-I_s-S_s

-2.80140879470855

GLOBAL_error=(1-(P_s+CH_in_m-(DRN_ext_m+dt_drn_m)-E_out_s-EVT_s-F_s-CH_out_m)/((S_to_m-S_from_m)+(S_f_s-S_i_s)))*100
print("GLOBAL Error:",GLOBAL_error, "%" )

GLOBAL Error: 0.5501517478422091 %

GLOBAL_error=P_s+CH_in_m-(DRN_ext_m+dt_drn_m)-E_out_s-EVT_s-F_s-CH_out_m-(S_to_m-S_from_m)-(S_f_s-S_i_s)
total_in=P_s+CH_in_m
percent_GLOBAL_error=GLOBAL_error/total_in
print("GLOBAL Error:",percent_GLOBAL_error, "%" )

GLOBAL Error: 0.002180232538103202 %

error_df=pd.DataFrame({"Coupling (DRN/INFLOW) % ": [error_drn_inflow], "Coupling (INF/RCH) %": [error_inf_rch], "SWMM %": [percent_SWMM_error], "MODFLOW %": [percent_MODFLOW_error], "Global %": [percent_GLOBAL_error]})

error_df

	Coupling (DRN/INFLOW) % 	Coupling (INF/RCH) % 	SWMM % 	MODFLOW % 	Global %
0 	-0.000041 	-0.000015 	-0.000173 	0.224647 	0.00218

Save Continuity Errors Data Frame

#error_df.to_excel("error_df_C2.xlsx")

SWMM Results
Time Series

    Subcatchments:

S_name="S4"
i=S_names_list.index(S_name)

f, (ax1, ax2, ax4, ax5) = plt.subplots(4,1, figsize=(10,13))
ax1.plot(S_SWMM_TS[i].Time,S_SWMM_TS[i]["Precipitation (mm/h)"], label=S_name,linestyle="--" )
ax1.set_ylabel("Precipitation (mm/h)")
ax1.legend()
ax2.plot(S_SWMM_TS[i].Time, S_SWMM_TS[i]["Infiltration (mm/h)"], label=S_name,linestyle="--"  )
ax2.set_ylabel("Infiltration (mm/h)")
ax2.legend()
# ax3.plot(S_SWMM_TS[i].Time,S_SWMM_TS[i]["Cumulative Infiltration (m3)"], label=S_name,linestyle="--" )
# ax3.set_ylabel("Cumulative Infiltration (m3)")
#ax3.legend()
ax4.plot(S_SWMM_TS[i].Time,S_SWMM_TS[i]["Runoff (m3/s)"], label=S_name,linestyle="--" )
ax4.legend()
ax4.set_ylabel("Runoff (m3/s)")
ax5.plot(S_SWMM_TS[i].Time,S_SWMM_TS[i]["Evaporation (mm/d)"], label=S_name,linestyle="--" )
ax5.legend()
ax5.set_ylabel("Evaporation (mm/d)")
f.suptitle("Subcatchment Time Series",fontsize=16)
plt.show()

    Storage Units:

WSU_name="W2"
i=WSU_names_list.index(WSU_name)

f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5,1, figsize=(10,12))
ax1.plot(WSU_SWMM_TS[i].Time,WSU_SWMM_TS[i]["Depth (m)"], label=WSU_name,linestyle="--" )
ax1.legend()
ax1.set_ylabel("Depth (m)")
ax2.plot(WSU_SWMM_TS[i].Time,WSU_SWMM_TS[i]["Total inflow (m3/s)"],label=WSU_name,linestyle="--" )
ax2.legend()
ax2.set_ylabel("Total inflow (m3/s)")
ax3.plot(WSU_SWMM_TS[i].Time,WSU_SWMM_TS[i]["Losses (m3/s)"],label=WSU_name,linestyle="--")
ax3.legend()
ax3.set_ylabel("Losses (m3/s)")
ax4.plot(WSU_SWMM_TS[i].index,WSU_SWMM_TS[i]["Lateral inflow (m3/s)"],label=WSU_name,linestyle="--"  )
ax4.legend()
ax4.set_ylabel("Latereal Inflow (m3/s)")
ax5.plot(WSU_SWMM_TS[i].Time,WSU_SWMM_TS[i]["Total outflow (m3/s)"],label=WSU_name,linestyle="--" )
ax5.legend()
ax5.set_ylabel("Total outflow (m3/s)")
f.suptitle("Wetland Time Series",fontsize=16)
plt.show()

FSU_name="SU1"
i=FSU_names_list.index(FSU_name)

f, (ax1, ax2, ax3, ax5) = plt.subplots(4,1, figsize=(10,12))
ax1.plot(FSU_SWMM_TS[i].Time,FSU_SWMM_TS[i]["Depth (m)"], label=FSU_name,linestyle="--" )
ax1.legend()
ax1.set_ylabel("Depth (m)")
ax2.plot(FSU_SWMM_TS[i].Time,FSU_SWMM_TS[i]["Total inflow (m3/s)"],label=FSU_name,linestyle="--" )
ax2.legend()
ax2.set_ylabel("Total inflow (m3/s)")
ax3.plot(FSU_SWMM_TS[i].Time,FSU_SWMM_TS[i]["Losses (m3/s)"],label=FSU_name,linestyle="--")
ax3.legend()
ax3.set_ylabel("Losses (m3/s)")
#ax4.plot(SU_SWMM_TS[i].index,SU_SWMM_TS[i]["Cumulative Exfiltration Loss (m3)"],label=SU_name,linestyle="--"  )
#ax4.legend()
#ax4.set_ylabel("Cumulative Exfiltration Loss (m3)")
ax5.plot(FSU_SWMM_TS[i].Time,FSU_SWMM_TS[i]["Total outflow (m3/s)"],label=FSU_name,linestyle="--" )
ax5.legend()
ax5.set_ylabel("Total outflow (m3/s)")
f.suptitle("Storage Unit Time Series",fontsize=16)
plt.show()

    Nodes:

J_name="TS21"
i=J_names_list.index(J_name)

f, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(10,12))
ax1.plot(J_SWMM_TS[i].Time,J_SWMM_TS[i]["Depth (m)"], label=J_name,linestyle="--" )
ax1.legend()
ax1.set_ylabel("Depth (m)")
ax2.plot(J_SWMM_TS[i].Time,J_SWMM_TS[i]["Flooding (m3/s)"],label=J_name,linestyle="--" )
ax2.legend()
ax2.set_ylabel("Flooding (m3/s)")
ax3.plot(J_SWMM_TS[i].Time,J_SWMM_TS[i]["Total inflow (m3/s)"],label=J_name,linestyle="--")
ax3.legend()
ax3.set_ylabel("Total inflow (m3/s)")
ax4.plot(J_SWMM_TS[i].Time,J_SWMM_TS[i]["Total outflow (m3/s)"],label=J_name,linestyle="--" )
ax4.legend()
ax4.set_ylabel("Total outflow (m3/s)")
f.suptitle("Junctions Time Series",fontsize=16)
plt.show()

    Weirs

V_name="O3"
i=V_names_list.index(V_name)

f, (ax1, ax2) = plt.subplots(2,1, figsize=(10,5))
ax1.plot(V_SWMM_TS[i].Time,V_SWMM_TS[i]["Flow (m3/s)"], label=V_name,linestyle="--" )
ax1.legend()
ax1.set_ylabel("Flow (m3/s)")
ax2.plot(V_SWMM_TS[i].Time,V_SWMM_TS[i]["Depth (m)"],label=V_name,linestyle="--" )
ax2.legend()
ax2.set_ylabel("Depth (m)")

Text(0, 0.5, 'Depth (m)')

MODFLOW Results
Spatial Results

    DRN_rate

from matplotlib.font_manager import FontProperties

plt.rcParams["font.family"] = "Times New Roman"
import matplotlib.colors as colors

dates_MODFLOW=pd.date_range(start="2017-03-01", end="2019-02-28", freq="D")

730

fig = plt.figure(figsize=(8, 8))
dia=(dates_MODFLOW == pd.Timestamp('2017-09-01')).argmax()
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax =TS_gdf[dia].plot(ax=ax, column='DRN_rate', legend=False, cmap='hot_r', vmin=0, vmax=600)
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
#ax.set_title("September 2017", fontsize=14)

pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
pmv.plot_grid(ax=ax, alpha=0.2);
plt.axis("off")
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams['mathtext.fontset'] = 'stix'
norm = colors.Normalize( vmin=0, vmax=600)
cbar = plt.cm.ScalarMappable(norm=norm, cmap='hot_r')
ax_cbar = fig.colorbar(cbar, ax=ax)
# add label for the colorbar
ax_cbar.set_label('Drain Rate $(\mathrm{m^3/day})$')

    RCH_rate (S)

fig = plt.figure(figsize=(8, 8))
dia=(dates_MODFLOW == pd.Timestamp('2017-09-01')).argmax()
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax =TS_gdf[dia].plot(ax=ax, column='RCH_S', legend=True, cmap='hot_r', vmin=0, vmax=1)
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
#ax.set_title("September 2017", fontsize=14)

pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
pmv.plot_grid(ax=ax, alpha=0.2);
plt.axis("off")
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams['mathtext.fontset'] = 'stix'
norm = colors.Normalize( vmin=0, vmax=1)
cbar = plt.cm.ScalarMappable(norm=norm, cmap='hot_r')
ax_cbar = fig.colorbar(cbar, ax=ax)
# add label for the colorbar
ax_cbar.set_label('Drain Rate $(\mathrm{m/day})$')

    Net Flux

dia=(dates_MODFLOW == pd.Timestamp('2018-09-01')).argmax()
net_diario=TS_gdf[dia].DRN_rate-(TS_gdf[dia].RCH_S.fillna(0)+TS_gdf[dia].RCH_WSU.fillna(0)+TS_gdf[dia].RCH_FSU.fillna(0))*x_resolution*y_resolution
gdf_diario=TS_gdf[dia]
gdf_diario["NET_FLUX"]=net_diario

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax =gdf_diario.plot(ax=ax, column='NET_FLUX', legend=False, cmap='seismic', vmax=400, vmin=-400)
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
#ax.set_title("September 2017", fontsize=14)

pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
pmv.plot_grid(ax=ax, alpha=0.2);
plt.axis("off")
plt.rc('font', size=14)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('legend', fontsize=14)    # legend fontsize
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams['mathtext.fontset'] = 'stix'
norm = colors.Normalize(vmax=400, vmin=-400)
cbar = plt.cm.ScalarMappable(norm=norm, cmap='seismic')
ax_cbar = fig.colorbar(cbar, ax=ax)
# add label for the colorbar
ax_cbar.set_label('Net Flux $(\mathrm{m^3/day})$')

Ground Water Heads

    LID Infiltration Criteria (1 m depth)

dia=(dates_MODFLOW == pd.Timestamp('2017-08-23')).argmax()
depth_dia=TS_gdf[dia].elevation-TS_gdf[dia].Altura
for i in range(len(depth_dia)):
    if type(TS_gdf[dia]["S"][i])!=str:
        depth_dia[i]=np.nan 
    if TS_gdf[dia].river[i] or TS_gdf[dia].lake[i]:
        depth_dia[i]=np.nan 
gdf_dia=TS_gdf[dia]
gdf_dia["depth"]=depth_dia
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

ax = gdf_dia.plot(ax=ax, column="depth",  cmap='RdYlGn', legend=False, vmax=3, vmin=-1)
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
#ax.set_title("August 2017", fontsize=14)
pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
pmv.plot_grid(ax=ax, alpha=0.2);
norm = colors.Normalize( vmax=3, vmin=-1)
cbar = plt.cm.ScalarMappable(norm=norm, cmap='RdYlGn')
ax_cbar = fig.colorbar(cbar, ax=ax)
ax_cbar.set_label('Depth Ground Water Table (m)')
plt.axis("off")

    Contour and Stream lines

dia=(dates_MODFLOW == pd.Timestamp('2017-09-1')).argmax()
x=np.linspace(min(np.array(MODFLOW_centroid['geometry'].x)),max(np.array(MODFLOW_centroid['geometry'].x)),100)
y=np.linspace(min(np.array(MODFLOW_centroid['geometry'].y)),max(np.array(MODFLOW_centroid['geometry'].y)),100)

z=list(TS_gdf[dia].Altura)
[x,y] = meshgrid(x,y)
z = np.array(z)
z = z.reshape((len(x), len(y)))
z=z[::-1]
fig, ax = plt.subplots(figsize=(8,8))
CS = ax.contour(x, y, z, 20)
ax.clabel(CS, inline=1, fontsize=10)
ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
dy, dx = np.gradient(-z) 
ax.streamplot(x, y, dx, dy, density=2,  color='0.9')

<matplotlib.streamplot.StreamplotSet at 0x14d44d92a08>

SAVE RESULTS

# def save_xls(list_dfs, xls_path):
#     with ExcelWriter(xls_path) as writer:
#         for n, df in enumerate(list_dfs):
#             df.to_excel(writer,'sheet%s' % n)
#         writer.save()

#save_xls(WSU_SWMM_TS,'Iteracion_mod_prec/WSU_SWMM_TS.xls')

#save_xls(FSU_SWMM_TS,'Iteracion_mod_prec/FSU_SWMM_TS.xls')

#save_xls(S_SWMM_TS,'Iteracion_mod_prec/S_SWMM_TS.xls')

#save_xls(J_SWMM_TS,'Iteracion_mod_prec/J_SWMM_TS.xls')

#save_xls(V_SWMM_TS,'Iteracion_mod_prec/V_SWMM_TS.xls')

#save_xls(C_SWMM_TS,'Iteracion_mod_prec/C_SWMM_TS.xls')

#save_xls(ZB_by_zone,'Iteracion_mod_prec/ZB_by_zone.xls')

#rch_drn_df=[S_rch_df, S_drn_df, WSU_rch_df, WSU_drn_df]

#save_xls(rch_drn_df,'Iteracion_mod_prec/rch_drn.xls')

Save monthly GDF

#dia=(dates_MODFLOW == pd.Timestamp('2018-06-01')).argmax()
#TS_gdf[dia].to_excel("GDF_results/2018_06_01.xlsx")

Determine Critic Storms

# rainfall_df=pd.DataFrame(index=S_SWMM_TS[0].Time)
# rainfall_df["mm"]=list(S_SWMM_TS[0]["Precipitation (mm/h)"])

# fig, ax = plt.subplots(figsize=(15,5))
# ax.axvline(x=pd.Timestamp('2017-08-23'), color='r')
# ax.axvline(x=pd.Timestamp('2018-02-15'), color='r')

# ax.plot(rainfall_df.index, rainfall_df.mm)
# plt.show()

Frecuency Analysis

# events_df=pd.DataFrame()
# events_df["Fecha"]=""
# events_df["Volumen [mm]"]=""
# events_df["Duracion [h]"]=""

# duracion=0
# mm=0
# i=0
# while i<(len(rainfall_df)):
#     if rainfall_df.mm[i]>0:
#         duracion+=1
#         mm+=rainfall_df.mm[i]
#         if rainfall_df.mm[i+1]==0:
#             events_df = events_df.append({'Volumen [mm]': mm , "Duracion [h]": duracion, "Fecha": rainfall_df.index[i]}, ignore_index=True)
#             duracion=0
#             mm=0
#             i+=1
#         else: 
#             i+=1
#     else:
#         i+=1

# df_events_sort=pd.DataFrame()
# df_events_sort=events_df.sort_values(by=["Volumen [mm]"], ascending=False)
# df_events_sort = df_events_sort.reset_index(drop=True)
# df_events_sort["Prob Excedencia"]=(df_events_sort.index+1)/(len(df_events_sort)+1) #corregir
# df_events_sort.head(20)

# dia=(dates_MODFLOW == pd.Timestamp('2017-08-23')).argmax()
# net_diario=TS_gdf[dia].DRN_rate-(TS_gdf[dia].RCH_S.fillna(0)+TS_gdf[dia].RCH_WSU.fillna(0)+TS_gdf[dia].RCH_FSU.fillna(0))*x_resolution*y_resolution
# gdf_diario=TS_gdf[dia]
# gdf_diario["NET_FLUX"]=net_diario

# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(1, 1, 1, aspect='equal')
# ax =gdf_diario.plot(ax=ax, column='NET_FLUX', legend=False, cmap='seismic', vmax=400, vmin=-400)
# ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
# ax.set_title('August 2017', fontsize=14)

# pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
# pmv.plot_grid(ax=ax, alpha=0.2);
# plt.axis("off")
# plt.rc('font', size=14)          # controls default text sizes
# plt.rc('axes', titlesize=14)     # fontsize of the axes title
# plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
# plt.rc('legend', fontsize=14)    # legend fontsize
# plt.rcParams["font.serif"] = ["Times New Roman"]
# plt.rcParams['mathtext.fontset'] = 'stix'
# norm = colors.Normalize(vmax=400, vmin=-400)
# cbar = plt.cm.ScalarMappable(norm=norm, cmap='seismic')
# ax_cbar = fig.colorbar(cbar, ax=ax)
# # add label for the colorbar
# ax_cbar.set_label('Net Flux $(\mathrm{m^3/day})$')
# #plt.savefig("Figures/Fig11a_August2017.svg")

# dia=(dates_MODFLOW == pd.Timestamp('2018-03-15')).argmax()
# TS_gdf[dia].to_excel("GDF_results/2018_03_15(lluvia_verano).xlsx")

# dia=(dates_MODFLOW == pd.Timestamp('2018-02-15')).argmax()
# net_diario=TS_gdf[dia].DRN_rate-(TS_gdf[dia].RCH_S.fillna(0)+TS_gdf[dia].RCH_WSU.fillna(0)+TS_gdf[dia].RCH_FSU.fillna(0))*x_resolution*y_resolution
# gdf_diario=TS_gdf[dia]
# gdf_diario["NET_FLUX"]=net_diario

# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(1, 1, 1, aspect='equal')
# ax =gdf_diario.plot(ax=ax, column='NET_FLUX', legend=False, cmap='seismic', vmax=400, vmin=-400)
# ax = SWMM_S_gdf.geometry.boundary.plot(color=None,edgecolor='k',linewidth = 0.5,ax=ax)
# ax.set_title('February 2018', fontsize=14)

# pmv = flopy.plot.PlotMapView(modelgrid=ml.modelgrid)
# pmv.plot_grid(ax=ax, alpha=0.2);
# plt.axis("off")
# plt.rc('font', size=14)          # controls default text sizes
# plt.rc('axes', titlesize=14)     # fontsize of the axes title
# plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
# plt.rc('legend', fontsize=14)    # legend fontsize
# plt.rcParams["font.serif"] = ["Times New Roman"]
# plt.rcParams['mathtext.fontset'] = 'stix'
# norm = colors.Normalize(vmax=400, vmin=-400)
# cbar = plt.cm.ScalarMappable(norm=norm, cmap='seismic')
# ax_cbar = fig.colorbar(cbar, ax=ax)
# # add label for the colorbar
# ax_cbar.set_label('Net Flux $(\mathrm{m^3/day})$')
# #plt.savefig("Figures/Fig11b_February2018.svg")

 

