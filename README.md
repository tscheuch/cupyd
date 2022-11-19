<p align="center">
  <img width="100" height="100" src="https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/apple/325/heart-with-arrow_1f498.png" alt="Cupyd">
</p>

<p align="center"><strong>Cupyd --</strong> a romantic integration scheme between SWMM & MODFLOW</p>

---

## Integration overview

The integration scheme done by **Cupyd** consists of three sequential steps:
1. A spatial integration between the [SWMM] and [MODFLOW] elements
2. A coupled model simulation with a spatio-temporal data exchange
3. Finally, a post processing with results analysis

### 1. Spatial integration between the SWMM and MODFLOW elements

**Input elements**
- SWMM model (swmm_model.inp) **or** a simulation object from [PySWMM]
- [Shapefiles](https://en.wikipedia.org/wiki/Shapefile) of
  SWMM subcatchments **and/or** storage units (S_polygon.shp **and/or** SU_polygon.shp)
  **and** junctions (J_points.shp)
- MODFLOW model **or** MODFLOW grid + DEM
  (only needed for spatial integration)
- **Optional:** List of groundwater zones and junctions’ association

**Output elements**
- [GeoDataFrame](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html)
  with MODFLOW cells’ associations (**elev**, **S**, **SU**, **DRN**, **drn_to**)
- A few plots

### 2. Coupled model simulation with a spatio-temporal data exchange

**Input elements**
- Simulation object from PySWMM
- [FloPy] model
- Time steps for reporting

**Output elements**
- Results from time series
- Simulation results from PySWMM
- [Zone budget](https://flopy.readthedocs.io/en/latest/source/flopy.utils.zonbud.html) results from FloPy

### 3. Post processing with results analysis

- Continuity analysis
- Results from SWMM
- Results from MODFLOW
- Integration results

## Library usage (work in progress)

For the setup and many other things, we have defined many utilities on the Makefile to facilitate your work.
None of them are quite complex, so you are welcome to go and check them or even run the commands separately.

1. Make sure that you have [Python](https://www.python.org) 3.9 (or higher) installed in your machine.

```sh
$ python --version
```

2. Make some [Poetry](https://python-poetry.org) by installing it using our Makefile target.

```sh
$ make poetry
```

3. Create a [virtual environment](https://docs.python.org/3/library/venv.html) to install the project dependencies.

```sh
$ make venv-with-dependencies
```

### Code health

In order to keep a healthy codebase, we use three different tools:
[black](https://github.com/psf/black),
[isort](https://github.com/PyCQA/isort) &
[mypy](https://github.com/python/mypy).

[modflow]:https://en.wikipedia.org/wiki/MODFLOW
[flopy]:https://github.com/modflowpy/flopy

[swmm]:https://en.wikipedia.org/wiki/Storm_Water_Management_Model
[pyswmm]:https://github.com/OpenWaterAnalytics/pyswmm
