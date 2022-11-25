from typing import Optional

import geopandas
import numpy

# from pyswmm import Simulation
from flopy.modflow import Modflow
from geopandas import GeoDataFrame


class CoupledModel:
    def __init__(
        self,
        modflow_model: Modflow,
        # pyswmm_model: Simulation,
        swmm_shp_file_path: str,
        storage_units_shp_file_path: Optional[str] = None,
        nodes_shp_file_path: Optional[str] = None,
    ) -> None:
        """
        This class intends to create the linkage (spatial integration) between a groundwater MODFLOW model and a surface SWMM model.
        It is a necessary and preliminary step for running a coupled SWMM-MODFLOW model.
        It allows the bidirectional and spatially distributed flux exchange between models, i.e. infiltration and exfiltration rates.

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
                (i.e., junctions, storage units, dividers or outfalls) that eventually can receive exfiltration
                as lateral inflow.

                Alternatively, if the relationship between MODFLOW cells and SWMM nodes is not built auotmatically,
                the file is a polygon shapefile that defines zones of the aquifer that drains to existing SWMM nodes.
                In both cases, one atribute for each geometry on the shapefile must be the name of the node in the SWMM model.
                By convention the name of the column should be "node".
                The shapefile must be previously built in a GIS software (QGIS, ArcGIS or similar).

        STEPS:
        1. Validation of MODFLOW model
            1.1 Validate `drn` package with only 1 stress period.
        2. Validation subcatchment-shape relationship is one-to-one.
            i.e. each existing subcatchment must have a polygon with the same name on the shapefile
        3. Make the spatial linktegration
        3.1. Create empty geodataframe
        3.2. Add MODFLOW info into empty dataframe
        3.3.
        """

        # self.validate_modflow_model()  # this validation should occur in the init
        # self.validate_subcatch_shp_relationship()  # this validation should occur in the init
        self.modflow_model = modflow_model
        self.swmm_shp_file_path = swmm_shp_file_path
        self.storage_units_shp_file_path = storage_units_shp_file_path
        self.nodes_shp_file_path = nodes_shp_file_path
        self.couple_models()

    def validate_modflow_model(self):
        """
        If the .dis file is missing, you will get an error since it's a required/minimal file for the coupling.
        If any additional file is missing, you will get a warning, since they are needed to run the simulation.
        You will not get an error, since it is still useful for coupling-related operation such as plotting.
        """
        pass

    def validate_subcatch_shp_relationship(self):
        pass

    def build_geodataframe(self) -> GeoDataFrame:
        """Function that builds and returns a `GeoDataFrame` with all the
        necessary columns for the spatial linkage between a `Modflow` model and a SWMM model.

        All the `GeoDataFrame` columns needed for the spacial linkage are:

        * MODFLOW-related columns:
            - x (int): Number of the x-axis were each cell lives in, on the `Modflow` instance grid.
            - y (int): Number of the y-axis were each cell lives in, on the `Modflow` instance grid.
            - geometry (Polygon): 'Polygon' representation of the cell.
            - elevation (float): Cell elevation of the first layer of the grid,
                the one that gets georeferenced with the SWMM model.
            - drn_elev (float): The elevation of the drain.
            - drn_cond (float): The hydraulic conductance of the interface between the aquifer and the drain.
            - ibound (int): It's the 1 if the cell is active; if not active, 0; if constant cells, -1.

        * SWMM-related columns:
            NOTE: un subcatchment va a infiltrar agua en múltiples celdas
                una celda recibe agua infiltrada desde un solo subcatchment
                una celda podría adicionalmente recibir agua desde un storage unit

            Infiltration exchange:
            - subcatchment (str): Name of SWMM subcatchment that infiltrates to the cell.
            - infiltration_storage_unit (str, None): Name of the SWMM storage unit that infiltrates to the cell.
            Exfiltration exchange:
            - node (str, None): Name of the node where the subcatchment exfiltrates. It can be a `junction`, a
                `storage unit`, `divider` or an `outfall`. These elements have unique names in between them, so it
                is safe to call them just `node`.

        TODO: Use a `temp` file for the generated modfow grid shapefile.
        TODO: Check how to avoid creating that file.

        See:
        https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/dis.html
        https://water.usgs.gov/ogw/modflow/MODFLOW-2005-Guide/drn.html

        Returns:
            GeoDataFrame: A `GeoDataFrame` instance with the described columns empty.
        """
        stress_period = 0
        top_layer_index = 0
        drn_elev = self.modflow_model.drn.stress_period_data.array["elev"][stress_period][
            top_layer_index
        ]

        drn_cond = self.modflow_model.drn.stress_period_data.array["cond"][stress_period][
            top_layer_index
        ]

        elevation = self.modflow_model.dis.top.array

        ibound = self.modflow_model.bas6.ibound[top_layer_index].array

        # TODO: Check what if some modflow model comes without projection
        # TODO: Dive deeper in `self.modflow_model.modelgrid.epsg`
        # The shape file adds the crs to the data frame
        self.modflow_model.modelgrid.write_shapefile("./temp_modflow.shp")
        self.geo_dataframe = geopandas.read_file("./temp_modflow.shp")
        self.geo_dataframe = self.geo_dataframe.drop(columns=["node"])
        self.geo_dataframe = self.geo_dataframe.rename(columns={"row": "x", "column": "y"})
        self.geo_dataframe["drn_elev"] = numpy.reshape(drn_elev, -1)
        self.geo_dataframe["drn_cond"] = numpy.reshape(drn_cond, -1)
        self.geo_dataframe["elevation"] = numpy.reshape(elevation, -1)
        self.geo_dataframe["ibound"] = numpy.reshape(ibound, -1)
        return self.geo_dataframe

    def couple_models(self):
        self.build_geodataframe()

        # Set the polygon centroids for the spatial joins
        self.geo_dataframe["centroids"] = self.geo_dataframe["geometry"].centroid
        self.geo_dataframe = self.geo_dataframe.set_geometry("centroids")

        # SPATIAL JOIN SUBCATCHMENTS

        swmm_geodf = geopandas.read_file(self.swmm_shp_file_path)

        joined_data = self.geo_dataframe.sjoin(swmm_geodf, how="inner", predicate="intersects")

        self.geo_dataframe["subcatchment"] = joined_data["S"]

        # SPATIAL JOIN STORAGE UNITS

        storage_unit_geodf = geopandas.read_file(self.storage_units_shp_file_path)

        joined_data = self.geo_dataframe.sjoin(
            storage_unit_geodf, how="inner", predicate="intersects"
        )

        self.geo_dataframe["infiltration_storage_unit"] = joined_data["stor_unit"]

        # SPATIAL JOIN NODES (only IF polygons)

        nodes_geodf = geopandas.read_file(self.nodes_shp_file_path)

        joined_data = self.geo_dataframe.sjoin(nodes_geodf, how="inner", predicate="intersects")

        self.geo_dataframe["node"] = joined_data["node"]

        self.geo_dataframe = self.geo_dataframe.set_geometry("geometry")
        return self.geo_dataframe
