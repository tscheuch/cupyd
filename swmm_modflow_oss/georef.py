import flopy
import geopandas as gpd

# from pyswmm import Simulation
from flopy.modflow import Modflow
from geopandas import GeoDataFrame
from shapely.geometry import Polygon


class CoupledModel:
    def __init__(
        self,
        modflow_model: Modflow,
        # pyswmm_model: Simulation,
        swmm_shp_file_path: str,
        storage_units_shp_file_path: str = None,
        nodes_shp_file_path: str = None,
    ) -> None:
        """
        This class intends to create the linkage (spatial integration) between a groundwater MODFLOW model and a surface SWMM model.
        It is a necessary and prelimary step for running a coupled SWMM-MODFLOW model. It allows the bidirectional and spatially distributed flux exchange between models (i.e., infiltration and exfiltration rates).

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

                TODO: falta averiguar cuáles son los archivos mínimos para correr una simulación flopy

            pyswmm_model (Simulation):
                <TODO: @tscheuch>

            subcatchments_swmm_shp_file_path (str):
                The filepath to the SWMM subcatchments shapefile.
                The SWMM subcatchments shapefile must be previously built in a GIS software (QGIS, ArcGIS or similar).
                It should be a polygon shapefile.
                One atribute for each geometry on the shapefile must be the name of the subcatchment in SWMM.
                By convention the name of the column should be "subcatch".

            storage_units_shp_file_path (str, optional):
                If the SWMM model has storage units with infiltration (seepage) capacity,
                to incorporated this rate as MODFLOW recharge, a storage unit shapefile must be given.
                The SWMM storage units shapefile must be previously built in a GIS software (QGIS, ArcGIS or similar).
                It represents the whole area associated to the respective water body.
                It should be a polygon shapefile.
                One atribute for each geometry on the shapefile must be the name of the storage unite in SWMM.
                By convention the name of the column should be "stor_unit".

            nodes_shp_file_path (str, optional):
                It is only necesary if the spatial linkage is in both directions, that is, if the drained
                water by the MODFLOW cells is incorporated as lateral inflow in SWMM nodes.

                If the relationship between MODFLOW cells and SWMM nodes is built auotmatically,
                the file is a point shapefile that represents the spatial location of the SWMM nodes
                (i.e., junctions, storage unitis, dividers or outfalls) that eventually can receive exfiltration
                as lateral inflow.

                Alternatively, if the relationship between MODFLOW cells and SWMM nodes is not built auotmatically,
                the file is a polygon shapefile that defines zones of the aquifer that drains to existing SWMM nodes.
                In both cases, one atribute for each geometry on the shapefile must be the name of the node in the SWMM model.
                By convention the name of the column should be "node".
                The shapefile must be previously built in a GIS software (QGIS, ArcGIS or similar).

        STEPS:
        1. Validation of modflow model
            1.1 Validate `drn` package with only 1 stress period.
        2. Validation subcatchment-shape relationship is one-to-one.
            i.e. each existing subcatchment must have a polygon with the same name on the shapefile
        3. Make the spatial linktegration
        3.1. Create empty geodataframe
        3.2. Add modflow info to empty dataframe
        3.3.
        """
        pass

        # self.validate_modflow_model()  # this validation should occur in the init
        # self.validate_subcatch_shp_relationship()  # this validation should occur in the init
        self.modflow_model = modflow_model

    def validate_modflow_model(self):
        """
        If the .dis file is missing, you will get an error since it's a required/minimal file for the coupling.
        If any additional file is missing, you will get a warning, since they are needed to run the simulation.
        You will not get an error, since it is still useful for coupling-related operation such as plotting.
        """
        pass

    def validate_subcatch_shp_relationship(self):
        pass

    @staticmethod
    def _create_empty_georeference_dataframe() -> GeoDataFrame:
        """Function that builds and returns an empty `GeoDataFrame` with all the
        necessary columns for the spatial linkage between a `Modflow` model and a SWMM model.

        All the `GeoDataFrame` columns needed for the spacial linkage are:

        * MODFLOW related columns:

            - x (Int): Number of the x-axis were each cell lives in, on the `Modflow` instance grid.
            - y (Int): Number of the y-axis were each cell lives in, on the `Modflow` instance grid.
            - geometry (Polygon): 'Polygon' representation of the cell.
            - elevation (Float): Cell elevation of the first layer of the grid,
                the one that gets georeferenced with the SWMM model.
            - drn_elev (Float): Is the elevation of the drain.
            - drn_cond (Float): Is the hydraulic conductance of the interface between the aquifer and the drain.

        * SWMM related columns:
            NOTE: un subcatchment va a infiltrar agua en múltiples celdas
                una celda recibe agua infiltrada desde un sólo subcatchment
                una celda podría adicionalmente recibir agua desde un storage unit

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
        geo_data_frame["x"] = ""
        geo_data_frame["y"] = ""
        geo_data_frame["geometry"] = ""
        geo_data_frame["elevation"] = ""
        geo_data_frame["drn_elev"] = ""
        geo_data_frame["drn_cond"] = ""

        # SWMM related columns
        geo_data_frame["subcatchment"] = ""
        geo_data_frame["infiltration_storage_unit"] = ""
        geo_data_frame["node"] = ""

        return geo_data_frame

    def _build_geometry_polygon(self, row, column):
        """Function to build a shapely geometry `Polygon` that represents the drain from
        the given row and column of the class `MODFLOW` model instance.

        For more information on the 'get_cell_vertices` of the modelgrid, see:
        https://flopy.readthedocs.io/en/3.3.5/source/flopy.discretization.structuredgrid.html?highlight=StructuredGrid#flopy.discretization.structuredgrid.StructuredGrid.get_cell_vertices

        Args:
            row (Int): The drain row index from the grid of the given model instance.
            column (Int): The drain column index from the grid of the given model instance.

        Returns:
            Polygon: A `Polygon` representing the drain.
        """
        # The function `get_cell_vertices` as the docs state:
        # "returns vertices for a single cell at row, column i, j."
        coords = self.modflow_model.modelgrid.get_cell_vertices(row, column)
        coords += [coords[0]]
        return Polygon(coords)

    def _build_modflow_related_entrance_data(self, row, column):
        """Function that returns an array containing the following information of a
        `MODFLOW` model instance top layer cell:

        - x (Int)
        - y (Int)
        - geometry (Polygon)
        - elevation
        - drn_elev (Float)
        - drn_cond (Float)

        See the `_create_empty_georeference_dataframe` docstrings for a better understanding on the variables.

        Args:
            row (Int): The drain row index from the grid of the given model instance.
            column (Int): The drain column index from the grid of the given model instance.

        Returns:
            List: A list containing all the described elements of a specific grid drain (cell)
                of the given `MODFLOW` model instance.
        """
        #  TODO: Find out what to do with multiple stress periods. Currently only using first stress period.
        stress_period = 0
        # TODO: @jricci1 Check the `stress_period_data` configuration. Done! Everything is on Telegram.
        # Need to decide where to write it.
        top_layer_index = 0
        geometry = self._build_geometry_polygon(row, column)
        elevation = self.modflow_model.dis.top.array[row][column]
        drn_elev = self.modflow_model.drn.stress_period_data.array["elev"][
            stress_period
        ][top_layer_index][row][column]
        drn_cond = self.modflow_model.drn.stress_period_data.array["cond"][
            stress_period
        ][top_layer_index][row][column]

        return [row, column, geometry, elevation, drn_elev, drn_cond]

    def _build_swmm_related_entrance_data(self):
        return ["", "", ""]

    def _build_data_frame_entrance(self, row, column):
        """Function that returns an array containing all the necesary information for the linkage
        of a specific `MODFLOW` cell with `SWMM` elements.

        See the `_create_empty_georeference_dataframe` docstrings for a better understanding on the variables.

        Args:
            modflow_model (Modflow): The `Modflow` model instance to get the data from.
            row (Int): The drain row index from the grid of the given model instance.
            column (Int): The drain column index from the grid of the given model instance.

        Returns:
            List: A list containing all the described elements of a specific grid drain (cell)
                of the given `MODFLOW` model instance.
        """
        modflow_data = self._build_modflow_related_entrance_data(row, column)
        swmm_data = self._build_swmm_related_entrance_data()
        return modflow_data + swmm_data

    def _add_modflow_info_to_dataframe(self, geodataframe) -> None:
        """Function that receives an empty `GeoDataFrame` to fill with its relevant
        `MODFLOW` model information.

        Information filled on the `GeoDataFrame`:
        - x (Int)
        - y (Int)
        - geometry (Polygon)
        - elevation (same as drn_elev)
        - drn_elev (Float)
        - drn_cond (Float)

        See the `_create_empty_georeference_dataframe` docstrings for a better understanding on the variables.

        ////////////
        IMPORTANT: WE NEED TO CONSIDER HOW IS THIS FUNCTION GOING TO WORK WITH MORE THAN ONE STRESS PERIOD (LAYER?).
        ////////////

        Args:
            geodataframe (GeoDataFrame): `DataFrame` to fill with information. It should come with
                the respective columns created.
        """
        model_grid_rows = self.modflow_model.modelgrid.nrow
        model_grid_cols = self.modflow_model.modelgrid.ncol

        for x in range(model_grid_rows):
            for y in range(model_grid_cols):
                cell = self._build_data_frame_entrance(x, y)
                geodataframe.loc[x * self.modflow_model.modelgrid.ncol + y] = cell


def empty_georeference_dataframe() -> GeoDataFrame:
    """Function that builds and returns an empty `GeoDataFrame` with all the
    necesary columns for the spatial linkage between a `Modflow` model and a SWMM model.

    All the `GeoDataFrame` columns needed for the georeference are:

    * MODFLOW related columns:

    ELIMINATE - node (Int): Unique identifier of each `Modflow` grid cell given by its row/column position.
        - row (Int): Number of the row were each cell lives in, on the `Modflow` instance grid.
        - column (Int): Number of the column were each cell lives in, on the `Modflow` instance grid.
        - geometry (Polygon): 'Polygon' representation of the cell.
        - elevation (Float): Cell elevation of the first layer of the grid,
            the one that gets georeferenced with the SWMM model.
        - drn_elev (Float): Is the elevation of the drain.
        - drn_cond (Float): Is the hydraulic conductance of the interface between the aquifer and the drain.

    * SWMM related columns:
        NOTE: un subcatchment va a infiltrar agua en múltiples celdas
              una celda recibe agua infiltrada desde un sólo subcatchment
              una celda podría adicionalmente recibir agua desde un storage unit

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
