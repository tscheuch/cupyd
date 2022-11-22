## Possible input-derived exceptions while using Cupyd

This document aims to define a non-exhaustive collection of anomalies
(i.e. errors and warnings) that could appear while using this library.
Please note that all of these exceptions **will originate from the input**.
Put simply, they should only appear if something is missing or ill-defined.

### Within the spatial integration module

#### Related to: spatial integration files

**ERROR: Missing MODFLOW files**
- **.dis** file is missing.

**WARNING: Missing MODFLOW files**
- **.drn** file is missing. Needed for a bidirectional flux exchange.
- **.bas** file is missing. Needed for a proper validation of the spatial integration.

**WARNING: Missing SWMM files**
- **.inp** file is missing. Needed for a proper validation of the GIS files.

**ERROR: Missing GIS files**
- **subcatchments_swmm_shp_file_path** is missing.

**WARNING: Missing GIS files**
- **storage_units_shp_file_path** is missing. Check **.inp** file.
- **nodes_shp_file_path** is missing. Needed for a bidirectional flux exchange.

**ERROR: Defective GIS files**
- **subcatchments_swmm.shp** must be a polygon shapefile.
- **storage_units.shp** must be a polygon shapefile.

**WARNING: Defective GIS files**
- **subcatchments_swmm.shp:** No column named "subcatch".
  First column is used by default.
- **storage_units.shp:** No column named "stor_unit".
  First column is used by default.
- **nodes.shp:** The file is a point shapefile.
  The relationship between MODFLOW cells and SWMM nodes is built automatically.
- **nodes.shp:** The file is a polygon shapefile.
  The relationship between MODFLOW cells and SWMM nodes is not built automatically.
- If **nodes.shp** is polygon: No column named "node" in **nodes.shp** file.
  First column is used by default.

**WARNING: Relation between GIS files and SWMM .inp file**
- Subcatchments in the **.inp** file == subcatchments_swmm.shp[subcatch]
- Storage units in the **.inp** file == storage_units.shp[stor_unit]
- Nodes in the **.inp** file == nodes.shp[node]
- Subcatchments areas in the **.inp** file == subcatchments_swmm.shp[subcatch] areas

#### Related to: spatial integration results

The following exceptions should be supported with spatial plots.

- MODFLOW spatial domain is smaller than the SWMM spatial domain.
  Part of the SWMM infiltraton will not be incorporated to MODFLOW as recharge. Continuity errors associated.
- There are SWMM subcatchments linked to inactive cells.
  SWMM infiltraton will not be incorporated to MODFLOW as recharge. Continuity errors associated.
- There are SWMM subcatchments linked to constant head cells.
  SWMM infiltraton will not be incorporated to MODFLOW as recharge. Continuity errors associated.
- There are cells with drain capacity that are not linked to any SWMM node. Continuity errors associated.
- The subcatchments areas in **subcatchments_swmm.shp** == sum cell areas associated to each subcatchment.
  It is recommended to use smaller cells.
- Inner cells are not associated to subcatchments.
- There are active cells outside the SWMM domain.

### Within the simulation module

- Inconsistent time steps
- SWMM flow units must be CMS
- Warning: SWMM simulation uses variable time steps
