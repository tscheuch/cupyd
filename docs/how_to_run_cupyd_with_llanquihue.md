# How to run Llanquihue with cupyd

This tutorial aims to explain how to use the cupyd library with the Llenquhue test that lives on the project.
It should cover everything so that a non technical person can jump into using it.

> This tutorial is build while the library resides on the [PyPi test repository](https://test.pypi.org/). It is expected be moved into PyPi at some point.

For this tutorial we do expect minimum terminal command knowledge. Basic command would be used.

Also, it assume that you have the following programs installed and in your PATH.

- Python >= 3.8
- Poetry >= 1.4


Files you need to run the example project.
- mfnwt --> MODFLOW executable program (For MAC/Linux. For windows this file should be called mfnwt.exe)
- Llanquihue
  - MODFLOW
    - LLANQUIHUE.bas
    - LLANQUIHUE.cbc
    - LLANQUIHUE.chx
    - LLANQUIHUE.dis
    - LLANQUIHUE.drn
    - LLANQUIHUE.hds
    - LLANQUIHUE.list
    - LLANQUIHUE.nam
    - LLANQUIHUE.nwt
    - LLANQUIHUE.oc
    - LLANQUIHUE.rch
    - LLANQUIHUE.upw
  - SWMM
    - GIS
      - SWMM_s --> shape related files
      - SWMM_SU --> shape related files
      - SWMM_nodes_zones --> shape related files
    - llanquihue_base.inp
    - llqh_horaria_16_19.dat

1. First move into the directory where we want to test the library.
   
```bash
cd your/project/directory
```

2. The we create the example project directory

```bash
mkdir example_project
```

3. Move into the directory

```bash
cd example_project
```

> Pro tip: Open VSCode in your project folder with: `code .`

4. Add the llanquihue, swmm and modflow specific files mentioned above

5. Create virtual environment

```bash
python -m venv .venv
```

6. Activate the virtual env

On Linux/Mac distributions
```bash
source .venv/bin/activate 
```

On Windows
```bash
.venv/Scripts/activate
```

This should create a `.env` folder in your current directory

7. Initialize the project with [poetry](https://python-poetry.org/docs/)

```
poetry init
```
This command will prompt you with a set of configurations for your project. For the sake of this tutorial, none are relevant. Feel free to just `enter` all the way to victory (it might require certain fields like the author).

> Poetry uses the `.venv` as the default virtual environment if present in the project directory.

--------------------------------------------------------------------
The following section is intended to be completed only while cupyd project lives in PyPi test repository.

Since the default repository to look for packages for poetry (as well as for pip) is PyPi, we need to change the default source repository so that it can find cupyd.

Add this piece of code to the `pyproject.toml`

```toml
[[tool.poetry.source]]
name = "test"
url = "https://test.pypi.org/simple/"
secondary = true
```
--------------------------------------------------------------------

8. Install cupyd

```bash
poetry add cupyd
```

9. Create the `llanquihu.py` file with

```python
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
```

8. Execute the simulation

```bash
python llanquihue.py
```
