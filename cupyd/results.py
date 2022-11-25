import pandas
from pyswmm import Links, Nodes, Simulation, Subcatchments


class AbstractSWMMObjectsTimeSeriesResult:
    """Abstract class to facilitate the creation and handling of
    different SWMM Objects time series result.

    The possible SWMM objects are:
    - Subcatchments
    - Storage units
    - Junctions
    - Links
    """

    def __init__(self) -> None:
        self._elements = {}

    @property
    def columns(self) -> list[str]:
        raise NotImplementedError("SWMMObjectsTimeSeriesColumns property must be defined")

    def build_empty_results_data_frame(self):
        return pandas.DataFrame(columns=self.columns)

    def add_empty_results(self, object_id):
        self._elements[object_id] = self.build_empty_results_data_frame()

    def __contains__(self, element_id):
        """
        Checks if Element ID exists.
        :return: ID Exists
        :rtype: bool
        """
        return element_id in self._elements

    def __getitem__(self, element_id):
        if self.__contains__(element_id):
            return self._elements[element_id]
        else:
            raise Exception("Element not in results")

    def __setitem__(self, element_id, name):
        self._elements[element_id] = name

    def __iter__(self):
        return self._elements.values()


class SubcatchmentsTimeSeriesResult(AbstractSWMMObjectsTimeSeriesResult):
    def __init__(self, simulation: Simulation) -> None:
        super().__init__()
        self._simulation = simulation
        self.add_simulation_subcatchments_empty_results(simulation)

    def add_simulation_subcatchments_empty_results(self, simulation: Simulation):
        for subcatchment in Subcatchments(simulation):
            self.add_empty_results(subcatchment.subcatchmentid)

    @property
    def columns(self):
        return [
            "Time",
            "Precipitation (mm/h)",
            "Evaporation (mm/d)",
            "Infiltration (mm/h)",
            "Runoff (m3/s)",
            "Runon (m3/s)",
            "Cumulative Infiltration (m3)",
            "Cumulative Evaporation (m3)",
        ]


class StorageUnitsTimeSeriesResult(AbstractSWMMObjectsTimeSeriesResult):
    def __init__(self, simulation: Simulation) -> None:
        super().__init__()
        self._simulation = simulation
        self.add_simulation_storage_units_empty_results(simulation)

    def add_simulation_storage_units_empty_results(self, simulation: Simulation):
        for node in Nodes(simulation):
            if node.is_storage():
                self.add_empty_results(node.nodeid)

    @property
    def columns(self):
        return [
            "Time",
            "Depth (m)",
            "Head (m)",
            "Flooding (m3/s)",
            "Lateral inflow (m3/s)",
            "Total inflow (m3/s)",
            "Total outflow (m3/s)",
            "Volume (m3)",
            "Losses (m3/s)",
            "Cumulative Exfiltration Loss (m3)",
            "Cumulative Evaporation Loss (m3)",
        ]


class JunctionsTimeSeriesResult(AbstractSWMMObjectsTimeSeriesResult):
    def __init__(self, simulation: Simulation) -> None:
        super().__init__()
        self._simulation = simulation
        self.add_simulation_junctions_empty_results(simulation)

    def add_simulation_junctions_empty_results(self, simulation: Simulation):
        for node in Nodes(simulation):
            if node.is_junction():
                self.add_empty_results(node.nodeid)

    @property
    def columns(self):
        return [
            "Time",
            "Depth (m)",
            "Head (m)",
            "Flooding (m3/s)",
            "Lateral inflow (m3/s)",
            "Total inflow (m3/s)",
            "Total outflow (m3/s)",
        ]


class LinksTimeSeriesResult(AbstractSWMMObjectsTimeSeriesResult):
    def __init__(self, simulation: Simulation) -> None:
        super().__init__()
        self._simulation = simulation
        self.add_simulation_links_empty_results(simulation)

    def add_simulation_links_empty_results(self, simulation: Simulation):
        for link in Links(simulation):
            if link.is_conduit():
                self.add_empty_results(link.linkid)

    @property
    def columns(self):
        return ["Time", "Depth (m)", "Head (m)", "Flooding (m3/s)"]
