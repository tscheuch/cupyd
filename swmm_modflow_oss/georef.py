import flopy
import geopandas as gpd
from flopy.modflow import Modflow
from geopandas import GeoDataFrame
from shapely.geometry import Polygon


class CoupledModel:
    def __init__(
        self,
        modflow_model: Modflow,
        swmm_shp_file_path: str,
        storage_units_shp_file_path: str = None,
        nodes_shp_file_path: str = None,
    ) -> None:
        """
        This class intends to create a relationship between a model subterraneo? in MODFLOW and
        a superficial one from SWMM.

        Args:
            modflow_model (Modflow):
                This is the `Modflow` model instance (which is based on flopy) to couple with a SWMM model.
                
                NOTE: Un modelo MODFLOW puede ser construido con distintas interfaces (i.e. no es necesario usar flopy)
                como por ejemplo Visual MODFLOW, porque finalmente lo relevante es generar este grupo de archivos (e.g. .dis, .bas, .drn)
                que en conjunto conforman el modelo MODFLOW. Lo que sí es relevante es poder importar todos estos archivos
                vía flopy para generar esta instancia que le entregaremos como argumento a esta clase.
                
                .dis --> Te da información acerca de la grilla/celdas: qué dimensiones (x,y,z) tiene, su elevación, cuántas capas, etc.
                         Además, te da la estructura inicial de la grilla.
                         El archivo .dis tiene forma matricial (row vs col) donde cada celda tiene info sobre la elevación.
                         Para sumar puntos en la trivia: .dis viene de discretization
                .drn --> Te dice qué celdas tienen capacidad de drenar, con qué capacidad (i.e. conductancia), a qué profundidad drenan
                .bas --> Viene de "basic". Te dice la actividad de las celdas: 1 -> activa, 0 -> inactiva, -1 -> celda de elevación constante
                
                If the .dis file is missing, you will get an error since it's a required/minimal file for the coupling.
                If any additional file is missing, you will get a warning, since they are needed to run the simulation.
                You will not get an error, since it is still useful for coupling-related operation such as plotting.
                

                TODO: falta averiguar cuáles son los archivos mínimos para correr una simulación flopy

            swmm_shp_file_path (str):
                The filepath to the SWMM `shp` file to load the externally generated file,
                that's usually created with a GIS-based software such as QGIS (open-source) or ARCGIS (proprietary).
                This `shp` file should be one made out of `Polygons`.
            
            storage_units_shp_file_path (str, optional):
                The filepath.
                This `shp` file should be one made out of `Polygons`.

            nodes_shp_file_path (str, optional):
                _description_.
        """
        pass


    def validate_modflow_model():
        """_summary_
        """

def empty_georeference_dataframe() -> GeoDataFrame:
    """Function that builds and returns an empty `GeoDataFrame` with all the
    necesary columns for the georeferenciation between a `Modflow` model and a SWMM model.

    All the `GeoDataFrame` columns needed for the georeference are:

    * MODFLOW related columns:

    ELIMINATE - node (Int): Unique identifier of each `Modflow` grid cell given by its row/column position.
        - row (Int): Number of the row were each cell lives in, on the `Modflow` instance grid.
        - column (Int): Number of the column were each cell lives in, on the `Modflow` instance grid.
        - geometry (Polygon): 'Polygon' representation of the cell.
        - elevation (Float): Cell elevation of the first layer of the grid,
            the one that gets georeferenced with the SWMM model. (Not necesary, seems to be the same as drn_elev)
        - drn_elev (Float): Is the elevation of the drain.
        - drn_cond (Float): Is the hydraulic conductance of the interface between the aquifer and the drain.

    * SWMM related columns:
        Infiltration exchange:
        - subcatchment (Str): Name of SWMM subcatchment that infiltrates to the cell.
        - infiltration_storage_unit (Str, None): Name of the SWMM storage unit that infiltrates to the cell.
        Exfiltration exchange:
        - node (Str, None): Name of the node where the subcatchment exfiltrate. It can be a `junction`, a
            `storage unit`, `divider` or an `outfall`. These elements have unique names in between them, so it
            is safe to call them just `node`.


    See:
    https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/dis.html
    https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/drn.html

    Returns:
        GeoDataFrame: A `GeoDataFrame` instance with the described columns empty.
    """
    geo_data_frame = GeoDataFrame()

    # MODFLOW related columns
    geo_data_frame["node"] = ""
    geo_data_frame["row"] = ""
    geo_data_frame["column"] = ""
    geo_data_frame["geometry"] = ""
    geo_data_frame["elevation"] = ""
    geo_data_frame["drn_elev"] = ""
    geo_data_frame["drn_cond"] = ""

    # SWMM related columns
    return geo_data_frame


def build_geometry_polygon(mf_model, row, column):
    """Function to build a shapely geometry `Polygon` that represents the drain from
    the given row and column of the given `MODFLOW` model instance.

    For more information on the 'get_cell_vertices` of the modelgrid, see:
    https://flopy.readthedocs.io/en/3.3.5/source/flopy.discretization.structuredgrid.html?highlight=StructuredGrid#flopy.discretization.structuredgrid.StructuredGrid.get_cell_vertices

    Args:
        mf_model (Modflow): The `Modflow` model instance to get the data from.
        row (Int): The drain row index from the grid of the given model instance.
        column (Int): The drain column index from the grid of the given model instance.

    Returns:
        Polygon: A `Polygon` representing the drain.
    """
    # The function `get_cell_vertices` as the docs state:
    # "returns vertices for a single cell at row, column i, j."
    coords = mf_model.modelgrid.get_cell_vertices(row, column)
    coords += [coords[0]]
    return Polygon(coords)


def build_data_frame_entrance(mf_model, stress_period, row, column):
    """Function that returns an array containing the following information of a
    `MODFLOW` model instance grid drain (cell):

    - node (Int)
    - row (Int)
    - column (Int)
    - geometry (Polygon)
    - elevation
    - drn_elev (Float)
    - drn_cond (Float)

    See the `empty_georeference_dataframe`docstrings for a better understanding on the variables.

    Args:
        mf_model (Modflow): The `Modflow` model instance to get the data from.
        stress_period (Int): The instance stress period from where to get the information.
        row (Int): The drain row index from the grid of the given model instance.
        column (Int): The drain column index from the grid of the given model instance.

    Returns:
        List: A list containing all the dscibed elements of a specific grid drain (cell)
            of the given `MODFLOW` model instance.
    """
    node = row * mf_model.modelgrid.ncol + column + 1
    geometry = build_geometry_polygon(mf_model, row, column)
    elevation = mf_model.dis.top.array[row][column]
    drn_elev = mf_model.drn.stress_period_data.array["elev"][stress_period][0][row][
        column
    ]  # drn_elev (Layer y stress period???)
    drn_cond = mf_model.drn.stress_period_data.array["cond"][stress_period][0][row][
        column
    ]  # drn_cond

    return [node, row, column, geometry, elevation, drn_elev, drn_cond]


def add_modflow_information_to_geodataframe(mf_model, geo_data_frame):
    """Function that receives a `MODFLOW` model instance and an empty `GeoDataFrame`
    to fill with its relevant information.

    Information filled on the `GeoDataFrame`:
    - node (Int)
    - row (Int)
    - column (Int)
    - geometry (Polygon)
    - elevation (same as drn_elev)
    - drn_elev (Float)
    - drn_cond (Float)

    See the `empty_georeference_dataframe`docstrings for a better understanding on the variables.

    ////////////
    IMPORTANT: WE NEED TO CONSIDER HOW IS THIS FUNCTION GOING TO WORK WITH MORE THAN ONE STRESS PERIOD (LAYER?).
    ////////////

    Args:
        mf_model (Modflow): The `Modflow` model instance to get the data from.
        geo_data_frame (GeoDataFrame): `DataFrame`to fill with information. It should come with
            the respective columns created.
    """
    model_grid_rows = mf_model.modelgrid.nrow
    model_grid_columns = mf_model.modelgrid.ncol

    for row in range(model_grid_rows):
        for col in range(model_grid_columns):
            entrance = build_data_frame_entrance(mf_model, 0, row, col)
            geo_data_frame.loc[row * mf_model.modelgrid.ncol + col] = entrance


def georeference_models(mf_model: Modflow, swmm_shp: str):
    """This function intends to couple the georeference between a MODFLOW
    model and the subcatchment elements of a SWMM model.

    Args:
        mf_model (Modflow): The `Modflow` model instance for the georeference.
        swmm_shp (str): The path to the SWMM `shp` file to load the handmade model
            subcatchemts georeference.

    Raises:
        TypeError: If the mf_model argument ain't a `Modflow` instance.

    Returns:
        GeoDataFrame: With all the columns filled with the georeferenciation
    """
    if not isinstance(mf_model, Modflow):
        raise TypeError(f"Attribute mf_model: {mf_model}, is not of type {Modflow}.")

    geo_data = empty_georeference_dataframe()
    add_modflow_information_to_geodataframe(mf_model, geo_data)

    # AGREGAR CRS --> importantísimo para la georeferenciación
    geo_data.crs = mf_model.modelgrid.proj4

    geo_data["centroids"] = geo_data["geometry"].centroid

    # AGREAGAR 2 columnas al gpd:
    # - subcuencas
    # - unidades de almacenamiento
    swmm_geodf = gpd.read_file(swmm_shp)

    geo_data = geo_data.set_geometry("centroids")

    joined_data = geo_data.sjoin(swmm_geodf, how="inner", predicate="intersects")

    joined_data = joined_data.set_geometry("geometry")
    return joined_data
