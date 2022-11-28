<p align="center">
  <img width="100" height="100" src="https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/320/apple/325/heart-with-arrow_1f498.png" alt="Cupyd">
</p>

<p align="center"><strong>Cupyd --</strong> a romantic integration scheme between SWMM & MODFLOW</p>

---

## Integration overview

The integration scheme done by **Cupyd** consists of three sequential steps:
1. A spatial integration between the [SWMM] and [MODFLOW] elements
2. A coupled model simulation with a spatio-temporal data exchange
3. Finally, a post-processing stage with results analysis

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
- A few [Matplotlib](https://matplotlib.org)-powered plots

### 2. Coupled model simulation with a spatio-temporal data exchange

**Input elements**
- Simulation object from PySWMM
- [FloPy] model
- Time steps for reporting

**Output elements**
- Results from time series
- Simulation results from PySWMM
- [Zone budget](https://flopy.readthedocs.io/en/latest/source/flopy.utils.zonbud.html) results from FloPy

### 3. Post-processing stage with results analysis

- Continuity analysis
- Results from SWMM
- Results from MODFLOW
- Integration results

## Library usage (work in progress)

We have defined a handful of custom Make targets in order to ease the local setup and development of Cupyd.
They’re quite simple, so you’re welcome to check them out and understand what’s going on behind the scenes.
Let’s use a few of them to run the example script.

1. Make sure that you have [Python](https://www.python.org) **3.9** (or higher) installed in your machine.

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

4. Run the [Llanquihue](https://en.wikipedia.org/wiki/Llanquihue_Lake)-based example script.

```sh
$ poetry run python llanquihue.py
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
