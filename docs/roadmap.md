# Roadmap del proyecto

The integration scheme is composed of 3 modules:
* Spatial Integration of SWMM and MODFLOW elements
* Coupled SWMM and MODFLOW Simulation: Loop for temporal and spatial data exchange
* Post processing and Results Analysis

## Spatial Integration of SWMM and MODFLOW elements

For this part, and on the first iteration, we want to build the `CoupledModel` class with the following capabilities:

- Initialize with a MODFLOW model instance, pyswmm Simulation instance , SWMM shp file path, storage unit shp file path (optional) and a nodes shp file path (optiona)
- Create a linkage beetween the models
- Allow plotting
