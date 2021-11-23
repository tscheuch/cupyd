# Integration scheme for coupling SWMM and MODFLOW
The integration scheme is composed of 3 modules:
* Spatial Integration of SWMM and MODFLOW elements
* Coupled SWMM and MODFLOW Simulation: Loop for temporal and spatial data exchange
* Post processing and Results Analysis
## Spatial Integration of SWMM and MODFLOW elements
### Inputs:
- SWMM model (swmm_model.inp) **or** Simulation from PySWMM
- SWMM subcatchments **or/and** storage units (S_polygon.shp **or/and** SU_polygon.shp) and junctions (J_points.shp) shapefiles
- MODLFOW model **or** MODFLOW grid + DEM (only for spatial integration)
- List of groundwater zones and junctiones association (optional) 
### Output: 
- gdf with MODFLOW cells associations (elev, S, SU, DRN, drn_to)
- plots
## Coupled SWMM and MODFLOW simulation: Loop for temporal and spatial data exchange 
### Inputs:
- Simulation from pySWMM
- MODFLOW flopy model (ml)
- reporting time steps 
### Outputs:
- Results Time Series
- Sim results (PySWMM)
- Zone Budjet (ZB) results (flopy)
## Post processing and Results analysis
- Continuity Analysis
- SWMM Results 
- MODFLOW Results
- Integration Results
