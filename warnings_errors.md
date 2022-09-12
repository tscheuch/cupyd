# Warnings and Errors associated to the Spatial Integration module
## Spatial Integration Files

### Error: Missing MODFLOW files
- Missing .dis file

### Warning: Missing MODFLOW files
- Missing .drn file. Rquired for a bidirectional flux exchange.
- Missing .bas file. Needed for suitable validation of the spatial integration.

### Warning: Missing SWMM files
- Missin .inp file. Needed for suitable validation of the GIS files.

### Error: Missing GIS files
- subcatchments_swmm_shp_file_path is missing.

### Warning: Missing GIS files
- storage_units_shp_file_path is missing (check .inp file)
- nodes_shp_file_path is missing. Needed for a bidirectional flux exchange.

### Warning: Wrong GIS files conceptualization
- subcatchments_swmm.shp: None column named "subcatch". First column is taken by default.
- storage_units.shp: None column named "stor_unit". First column is taken by default.
- nodes.shp: The file is a point shapefile. The relationship between MODFLOW cells and SWMM nodes is built auotmatically.
- nodes.shp: The fie is a polygon shapefile. The relationship between MODFLOW cells and SWMM nodes is not built auotmatically.
- (if nodes.shp is polygon): None column named "node" in nodes.shp file. First column is taken by default.

### Error: Wrong GIS files conceptualization
- subcatchments_swmm.shp must be a polygon shapefile.
- storage_units.shp must be a polygon shapefile.
### Warning: Relation between GIS files and SWMM .inp file
- subcatchments in the .inp file == subcatchments_swmm.shp["subcath"]
- storage uniste in the .inp file== storage_units.shp [stor_unit]
- nodes in the .inp file == nodes.shp[node]
- subcatchments_swmm.shp["subcath"] areas == subcatchments areas in the .inp file


## Spatial Integration Results
(This warnings and errors should be supported with spatial plots)

*Warning:
- MODFLOW spatial domain is smaller tha SWMM spatial domain. Part of the SWMM infiltraton will not be incorporated to MODFLOW as recharge. Continuity errors associated. 
- There are SWMM subcathments linked to inactive cells. SWMM infiltraton will not be incorporated to MODFLOW as recharge. Continuity errors associated. 
- There are SWMM subcathments linked to constant head cells. SWMM infiltraton will not be incorporated to MODFLOW as recharge. Continuity errors associated. 
- There are cells with drain capacity that are not linked to any SWMM node. Continuity errors associated. 
- The subcathments areas in subcatchments_swmm.shp == sum cell areas asociated to each subcathcment. It is recommended to use smaller cells.
- Inner cells are not associated to subcatchments.
- There are active cells outside the SWMM domain. 

# Warnings and Errors associated to the Simulation module
- Consistent times steps
- SWMM flow units must be CMS
- Warning: SWMM simulation uses variable time steps