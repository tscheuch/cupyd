<p align="center">
  <img width="100" height="100" src="https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/apple/325/heart-with-arrow_1f498.png" alt='cupyd'>
</p>

<p align="center"><strong>cupyd --</strong> a romantic integration scheme between SWMM & MODFLOW</p>

---

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

# Library usage (WIP)

For the setup and many other things, we have defined many utilities on the Makefile to facilitate your work.
None of them are quite complex, so you are welcome to go and check them or even run the commands separately.

To facilitate the instructions, we use them all around.

1. Ensure that you have Python 3.8 (or higher) installed.

2. Make [Poetry](https://python-poetry.org/) by installing it with our custom Makefile target.
```sh
$ make poetry
```

3. Create a [virtual environment](https://docs.python.org/3/library/venv.html) to install the project dependencies.
```sh
$ make venv-with-dependencies
```
