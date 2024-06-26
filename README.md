<p align="center">
  <img width="100" height="100" src="https://em-content.zobj.net/thumbs/240/apple/354/heart-with-arrow_1f498.png" alt="Cupyd">
</p>

<p align="center"><strong>Cupyd --</strong> a romantic integration scheme between SWMM & MODFLOW</p>

---

## Integration overview

The integration scheme done by **Cupyd** consists of two sequential steps:
1. A spatial integration between the [SWMM] and [MODFLOW] elements
2. A coupled model simulation with a spatio-temporal data exchange


### 1. Spatial integration between the SWMM and MODFLOW elements

Code file: georef.py

**Input elements**
- SWMM model (swmm_model.inp) **or** a simulation object from [PySWMM]
- [Shapefiles]  mandatory files (SHP, SHX and DBF) (https://en.wikipedia.org/wiki/Shapefile) of:
  SWMM subcatchments **and/or** SWMM storage units
  **and**, if a bidirecational exchange is needed, SWMM nodes. Note that a distinction has been made between nodes and storage units, as the latter are the only ones with drainage capacity. In the nodes shapefiles storage units can be incorporated.

All shapefiles must be of polygon type. In the case of subcatchments and storage units, they represent the surface area they occupy in space. For junctions, they represent the area of the aquifer that will eventually drain water to that point.

To prevent continuity errors in the flow transfer, it's important to avoid overlapping polygons and empty spaces. Additionally, if you want bidirectional exchange, make sure that all MODFLOW cells with drain capacity intersect with a polygon from the junction shapefile.

One atribute for each geometry on the shapefile must be the name of the element in SWMM .inp file.
By convention the name of the column should be "subcatch", "stor_unit" and "node".


- MODFLOW model **or** MODFLOW grid shapefile + Digital Elevation Model (DEM)
  (only needed for spatial integration)


**Output elements**

- [GeoDataFrame](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html)
  with MODFLOW cells’ associations (**x**,**y**,**geometry**, **elev**, **drn_elev**, **drn_cond** , **subcatch**, **stor_unit**, **drn_to_node**, **drn_rate**), where:

  **x**: Number of the x-axis were each cell lives in, on the `Modflow` instance grid.  
  **y**: Number of the y-axis were each cell lives in, on the `Modflow` instance grid.  
  **geometry**: 'Polygon' representation of the cell.  
  **elev**: Cell elevation of the first layer of the grid, the one that gets georeferenced with the SWMM model.  
  **drn_elev**: The elevation of the drain.  
  **drn_cond**: The hydraulic conductance of the interface between the aquifer and the drain.  

  **subcatch**: Name of the SWMM subcatchment that infiltrates to the cell.  
  **stor_unit**: Name of the SWMM storage unit that infiltrates to the cell.  
  **drn_to_node**: Name of the node where the cell eventually exfiltrates. It can be a `junction`, a `storage unit`, `divider` or an `outfall`.  

  Note that a subcatchment and a storage unit have the ability to recharge multiple cells. Moreover, a single cell can receive infiltration from both a subcatchment and a storage unit simultaneously. However, in the opposite direction, a cell can only drain to one specific node.

- A few [Matplotlib](https://matplotlib.org)-powered plots where you can visually verify if the spatial integration was correctly constructed.

### 2. Coupled model simulation with a spatio-temporal data exchange

**Input elements**
- Simulation from PySWMM
- [FloPy] model
- Spatial integration Geodataframe  

**Output elements**
- Geodataframe with exchanged flux can be saved/ploted in every time step
- Simulation results from PySWMM
- [Zone budget](https://flopy.readthedocs.io/en/latest/source/flopy.utils.zonbud.html) results from FloPy (developing)


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
