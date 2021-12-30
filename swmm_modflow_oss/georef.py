import flopy

import geopandas as gpd

from shapely.geometry import Polygon


def georeference_models(mf_model=None, swmm_shp=None):
    if not isinstance(mf_model, flopy.modflow.Modflow):
        raise TypeError(f"Attribute mf_model: {mf_model}, is not of type {flopy.modflow.Modflow}.")
        
    geo_data = gpd.GeoDataFrame()
    geo_data["node"] = ""
    geo_data["row"] = ""
    geo_data["column"] = ""
    geo_data["geometry"] = ""
    geo_data["elevation"] = ""
    geo_data["drn_elev"] = ""
    geo_data["drn_cond"] = ""
    node = 0
    rows = mf_model.modelgrid.nrow
    cols = mf_model.modelgrid.ncol
    for row in range(rows):
        for col in range(cols):
            coords = mf_model.modelgrid.get_cell_vertices(row, col,0)
            coords += [coords[0]]
            polygon = Polygon(coords)
            node += 1
            entrance = [
                row * mf_model.modelgrid.ncol + col + 1, # node
                row + 1, # row
                col + 1, # col
                polygon, # geometry
                mf_model.dis.top.array[row][col], # elevation
                mf_model.drn.stress_period_data.array["elev"][0][0][row][col], # drn_elev (Layer y stress period???)
                mf_model.drn.stress_period_data.array["cond"][0][0][row][col] # drn_cond
            ]
            geo_data.loc[row * mf_model.modelgrid.ncol + col] = entrance

    # AGREGAR CRS --> importantísimo para la georeferenciación
    geo_data.crs = mf_model.modelgrid.proj4

    geo_data["centroids"] = geo_data["geometry"].centroid

    # AGREAGAR 2 columnas al gpd:
    # - subcuencas
    # - unidades de almacenamiento
    swmm_geodf = gpd.read_file(swmm_shp)

    geo_data = geo_data.set_geometry('centroids')

    joined_data = geo_data.sjoin(swmm_geodf, how="inner", predicate='intersects')

    joined_data = joined_data.set_geometry('geometry')
    return joined_data
