"""
This module provides classes and methods for standardizing and processing GeoDataFrames.

It includes functionality for transforming input GeoDataFrames or file paths into standardized formats,
removing Z-coordinates, handling invalid geometries, and processing polygons to remove "geni" (holes).
It supports both GeoPandas frameworks.

Author: Sergio A. Pelayo Escalera
Email: sergioapelayoe@gmail.com
Created: 2025-04-09
Last-modified: 2025-04-09
Version: 0.1.0
"""

__version__ = "0.1.0"
__author__ = "Sergio A. Pelayo Escalera (sergioapelayoe@gmail.com)"
__collabs__ = [""]
__created__ = "2025-04-09"
__last_updated__ = "2025-04-09"

import os
import numpy as np
from typing import Union, List, Tuple, Optional, Any
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon, LinearRing, GeometryCollection


class GeoDataFrameAdapter:
    """
    Adapter class to provide a unified interface for GeoPandas GeoDataFrames.

    This adapter wraps different GeoDataFrame implementations to provide consistent access
    to their functionality, regardless of the underlying framework used.
    """

    def __init__(self,
                 geodataframe: Any):
        """
        Initialize the GeoDataFrameAdapter with a GeoDataFrame.

        Args:
            geodataframe: A GeoPandas GeoDataFrame GeoDataFrame

        Raises:
            TypeError: If the provided object is not a supported GeoDataFrame type
        """
        self.geodataframe = geodataframe
        if isinstance(geodataframe, gpd.GeoDataFrame):
            self.framework = "geopandas"
        else:
            raise TypeError(
                "Input must be a GeoPandas or GeoPolars GeoDataFrame")

    def to_crs(self,
               crs: str,
               inplace: bool = False) -> Any:
        """
        Transform the GeoDataFrame to a new coordinate reference system.

        Args:
            crs: The coordinate reference system to transform to
            inplace: Whether to modify the GeoDataFrame in-place

        Returns:
            GeoDataFrame: The transformed GeoDataFrame
        """
        if self.framework == "geopandas":
            if inplace:
                self.geodataframe.to_crs(crs, inplace=True)
                return self.geodataframe
            else:
                return self.geodataframe.to_crs(crs)

    @property
    def crs(self) -> Optional[str]:
        """
        Get the Coordinate Reference System (CRS) of the GeoDataFrame.

        Returns:
            Optional[str]: The CRS of the GeoDataFrame, or None if not set.
        """
        return self.geodataframe.crs

    @crs.setter
    def crs(self,
            value: str) -> None:
        """
        Set the Coordinate Reference System (CRS) of the GeoDataFrame.

        Args:
            value (str): The CRS to set, typically in the format "EPSG:XXXX".
        """
        self.geodataframe.crs = value

    def copy(self) -> Any:
        """
        Create a deep copy of the GeoDataFrame.

        Returns:
            Any: A new instance with a copy of the underlying GeoDataFrame.
        """
        return self.geodataframe.copy()

    @property
    def columns(self) -> List[str]:
        """
        Get the column names of the GeoDataFrame.

        Returns:
            List[str]: A list containing the column names of the GeoDataFrame.
        """
        return list(self.geodataframe.columns)

    def __getitem__(self,
                    key):
        """
        Get items from the GeoDataFrame using indexing or column names.

        Args:
            key: The index, column name, or slice to retrieve.

        Returns:
            The selected data from the GeoDataFrame.
        """
        return self.geodataframe[key]

    def iterrows(self):
        """
        Iterate over rows of the GeoDataFrame as (index, Series) pairs.

        Returns:
            Iterator: An iterator yielding pairs of (index, Series) for each row.
        """
        return self.geodataframe.iterrows()

    def apply(self,
              func,
              axis=0) -> Any:
        """
        Apply a function along an axis of the GeoDataFrame.

        Args:
            func: Function to apply to each column or row.
            axis (int, default 0): Axis along which the function is applied.
                0 or 'index': apply function to each column.
                1 or 'columns': apply function to each row.

        Returns:
            The result of applying the function to the GeoDataFrame.
        """
        return self.geodataframe.apply(func, axis=axis)

    def to_native(self) -> Any:
        """
        Return the native GeoDataFrame object.

        This method provides access to the underlying GeoPandas GeoDataFrame
        for operations that require the native object.

        Returns:
            Any: The underlying GeoPandas GeoDataFrame object.
        """
        return self.geodataframe

    @classmethod
    def from_file(cls,
                  path: str) -> 'GeoDataFrameAdapter':
        """
        Create a GeoDataFrameAdapter from a file.

        Args:
            path: Path to the file

        Returns:
            GeoDataFrameAdapter: A new adapter instance
        """
        if not os.path.exists(path):
            raise ValueError(f"Path {path} does not exist.")
        gdf = gpd.read_file(path)
        return cls(gdf)


class _GeniRemover:
    """
    A class to remove "geni" (holes) from a Polygon geometry.

    This class processes a Polygon geometry by identifying and merging its holes
    (interior rings) with the exterior ring, ensuring the resulting geometry has no holes.

    Attributes:
        geom (Polygon): The input Polygon geometry to process.
    """

    def __init__(self,
                 geom: Polygon):
        """
        Initialize the _GeniRemover instance.

        Args:
            geom (Polygon): The input Polygon geometry to process.
        """
        self.geom = geom

    def _calculate_distance(self,
                            point1: Tuple[float, float],
                            point2: Tuple[float, float]) -> float:
        """
        Calculate the Euclidean distance between two points.

        Args:
            point1 (Tuple[float, float]): The first point (x, y).
            point2 (Tuple[float, float]): The second point (x, y).

        Returns:
            float: The Euclidean distance between the two points.
        """
        return np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

    def _initialize_nearest_search(self) -> Tuple[None, float]:
        """
        Initialize the nearest search with default values.

        Returns:
            Tuple[None, float]: A tuple containing None and infinity as the initial values.
        """
        return None, float('inf')

    def _update_nearest_pair(self,
                             ext_point: Tuple[float, float],
                             hole_point: Tuple[float, float],
                             nearest_pair: Tuple[Tuple[float, float], Tuple[float, float]],
                             min_dist: float) -> Tuple[Tuple[Tuple[float, float], Tuple[float, float]], float]:
        """
        Update the nearest pair of points between the exterior and a hole.

        Args:
            ext_point (Tuple[float, float]): A point from the exterior ring.
            hole_point (Tuple[float, float]): A point from the hole.
            nearest_pair (Tuple[Tuple[float, float], Tuple[float, float]]): The current nearest pair of points.
            min_dist (float): The current minimum distance.

        Returns:
            Tuple[Tuple[Tuple[float, float], Tuple[float, float]], float]: The updated nearest pair and distance.
        """
        dist = self._calculate_distance(ext_point, hole_point)
        if dist < min_dist:
            return (ext_point, hole_point), dist
        return nearest_pair, min_dist

    def _find_nearest_points(self,
                             exterior_coords: List[Tuple[float, float]],
                             hole_coords: List[Tuple[float, float]]) -> Tuple[Tuple[float, float], Tuple[float, float], float]:
        """
        Find the nearest pair of points between the exterior and a hole.

        Args:
            exterior_coords (List[Tuple[float, float]]): Coordinates of the exterior ring.
            hole_coords (List[Tuple[float, float]]): Coordinates of the hole.

        Returns:
            Tuple[Tuple[float, float], Tuple[float, float], float]: The nearest pair of points and their distance.
        """
        nearest_pair, min_dist = self._initialize_nearest_search()
        for ext_point in exterior_coords:
            for hole_point in hole_coords:
                nearest_pair, min_dist = self._update_nearest_pair(ext_point,
                                                                   hole_point,
                                                                   nearest_pair,
                                                                   min_dist)
        return nearest_pair[0], nearest_pair[1], min_dist

    def _process_hole(self,
                      curr_ext: List[Tuple[float, float]],
                      hole: List[Tuple[float, float]],
                      ext_point: Tuple[float, float],
                      hole_point: Tuple[float, float]) -> List[Tuple[float, float]]:
        """
        Merge a hole into the exterior ring by creating a cut.

        This method creates a "cut" from the exterior ring to the hole and back,
        effectively removing the hole by incorporating it into the exterior boundary.
        The hole is traversed in reverse (clockwise) to ensure proper polygon topology.

        Args:
            curr_ext (List[Tuple[float, float]]): The current exterior ring coordinates.
            hole (List[Tuple[float, float]]): The hole coordinates to be merged.
            ext_point (Tuple[float, float]): The nearest point on the exterior ring to the hole.
            hole_point (Tuple[float, float]): The nearest point on the hole to the exterior ring.

        Returns:
            List[Tuple[float, float]]: The updated exterior ring coordinates with the hole merged.
        """
        insert_idx = curr_ext.index(ext_point)
        hole_point_idx = hole.index(hole_point)
        # The key change: traverse the hole in REVERSE direction (clockwise)
        # This ensures we're "cutting out" the hole rather than adding it
        ordered_hole = ([hole_point] +
                        # Reverse from hole_point to start
                        hole[hole_point_idx-1::-1] +
                        hole[:hole_point_idx-1:-1])   # Reverse from end to hole_point
        # Create a cut from exterior to hole and back
        new_sequence = [ext_point] + ordered_hole + [ext_point]
        # Insert the cut sequence into the exterior
        return curr_ext[:insert_idx] + new_sequence + curr_ext[insert_idx+1:]

    def _find_nearest_hole(self,
                           curr_ext: List[Tuple[float, float]],
                           holes: List[List[Tuple[float, float]]]) -> Tuple[int, Tuple[float, float], Tuple[float, float]]:
        """
        Find the nearest hole to the exterior ring.

        This method identifies the hole that is closest to the exterior ring by calculating
        the minimum distance between points on the exterior and points on each hole.

        Args:
            curr_ext (List[Tuple[float, float]]): The current exterior ring coordinates.
            holes (List[List[Tuple[float, float]]]): A list of hole coordinates.

        Returns:
            Tuple[int, Tuple[float, float], Tuple[float, float]]:
                - The index of the nearest hole in the `holes` list.
                - The nearest point on the exterior ring.
                - The nearest point on the hole.
        """
        min_dist = float('inf')
        nrst_hole_index = None
        nrst_ext_point = None
        nrst_hole_point = None
        for idx, hole in enumerate(holes):
            ext_point, hole_point, dist = self._find_nearest_points(
                curr_ext, hole)
            if dist < min_dist:
                min_dist = dist
                nrst_hole_index = idx
                nrst_ext_point = ext_point
                nrst_hole_point = hole_point
        return nrst_hole_index, nrst_ext_point, nrst_hole_point

    def _merge_holes_with_exterior(self,
                                   holes: List[List[Tuple[float, float]]],
                                   curr_ext: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
        """
        Merge all holes into the exterior ring.

        Args:
            holes (List[List[Tuple[float, float]]]): A list of hole coordinates.
            curr_ext (List[Tuple[float, float]]): The current exterior ring coordinates.

        Returns:
            List[Tuple[float, float]]: The updated exterior ring coordinates.
        """
        while holes:
            nrst_hole_index, nrst_ext_point, nrst_hole_point = self._find_nearest_hole(
                curr_ext, holes)
            curr_ext = self._process_hole(curr_ext,
                                          holes[nrst_hole_index],
                                          nrst_ext_point,
                                          nrst_hole_point)
            holes.pop(nrst_hole_index)
        return curr_ext

    def _close_ring(self,
                    coords: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
        """
        Ensure the ring is closed by appending the first coordinate to the end.

        Args:
            coords (List[Tuple[float, float]]): The ring coordinates.

        Returns:
            List[Tuple[float, float]]: The closed ring coordinates.
        """
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        return coords

    def trnsf_pol_all_geni(self) -> Polygon:
        """
        Transform the Polygon by removing all holes.

        Returns:
            Polygon: The transformed Polygon without holes.
        """
        if not self.geom.interiors:
            return self.geom
        holes = [list(interior.coords) for interior in self.geom.interiors]
        curr_ext = list(self.geom.exterior.coords)
        curr_ext = self._merge_holes_with_exterior(holes, curr_ext)
        curr_ext = self._close_ring(curr_ext)
        return Polygon(curr_ext)


class _InputTransformer:
    """
    A class to transform and standardize input GeoDataFrames or file paths.

    This class handles input validation, CRS standardization, and geometry processing
    to remove Z-coordinates and "geni" (holes) from polygons.

    Attributes:
        gdf_input (Union[gpd.GeoDataFrame, str]): The input GeoDataFrame or file path.
        gdf (gpd.GeoDataFrame): The standardized GeoDataFrame.
    """

    def __init__(self,
                 gdf_input: Union[gpd.GeoDataFrame, str, Any],
                 crs: str = "EPSG:4326",
                 remove_geni: bool = True) -> None:
        """
        Initialize the _InputTransformer instance.

        Args:
            gdf_input (Union[gpd.GeoDataFrame, str, Any]): The input Geopandas GeoDataFrame, or file path.
            crs (str, optional): The coordinate reference system to use. Defaults to "EPSG:4326".

        Raises:
            ValueError: If the input is neither a valid file path nor a supported GeoDataFrame.
        """
        self.gdf_input = gdf_input
        self.target_crs = crs
        self.remove_geni = remove_geni
        self.standard_gdf = self.apply_z_coord_and_geni_removal(
            self.set_standard_crs())

    def _is_path(self) -> bool:
        """
        Check if the input is a file path.

        Returns:
            bool: True if the input is a file path, False otherwise.
        """
        return isinstance(self.gdf_input, str)

    def _is_geodataframe(self) -> bool:
        """
        Check if the input is a supported GeoDataFrame.

        Returns:
            bool: True if the input is a supported GeoDataFrame, False otherwise.
        """
        is_gpd = isinstance(self.gdf_input, gpd.GeoDataFrame)
        return is_gpd

    def _path_validation(self) -> None:
        """
        Validate the file path.

        Raises:
            ValueError: If the file path does not exist.
        """
        if not os.path.exists(self.gdf_input):
            raise ValueError(f"Path {self.gdf_input} does not exist.")

    def _get_gdf(self) -> gpd.GeoDataFrame:
        """
        Retrieve the GeoDataFrame from the input.

        Returns:
            gpd.GeoDataFrame: The GeoDataFrame.

        Raises:
            ValueError: If the input is neither a valid file path nor a supported GeoDataFrame.
        """
        if self._is_path():
            self._path_validation()
            return gpd.read_file(self.gdf_input)
        elif self._is_geodataframe():
            if isinstance(self.gdf_input, gpd.GeoDataFrame):
                return self.gdf_input
        else:
            raise ValueError(
                "Input should be a path or a GeoPandas GeoDataFrame.")

    def set_standard_crs(self) -> gpd.GeoDataFrame:
        """
        Standardize the CRS of the GeoDataFrame to the target CRS.

        Returns:
            gpd.GeoDataFrame: The GeoDataFrame with a standardized CRS.
        """
        gdf = self._get_gdf()
        if gdf.crs is None:
            gdf.crs = self.target_crs
        gdf.to_crs(self.target_crs, inplace=True)
        return gdf

    def _if_linear_ring_to_polygon(self,
                                   geom) -> Polygon:
        """
        Convert a LinearRing to a Polygon if applicable.

        Args:
            geom: The geometry to check.

        Returns:
            Polygon: The converted Polygon or the original geometry.
        """
        if isinstance(geom, LinearRing):
            return Polygon(geom)
        return geom

    def _if_collection_of_linear_rings_to_multipolygon(self,
                                                       geom) -> MultiPolygon:
        """
        Convert a GeometryCollection of LinearRings to a MultiPolygon.

        Args:
            geom: The geometry to check and potentially convert.

        Returns:
            MultiPolygon: The converted MultiPolygon if the input is a GeometryCollection
            of LinearRings, otherwise the original geometry is returned unchanged.

        Raises:
            Exception: If the GeometryCollection is empty or not composed entirely of
                valid LinearRings.
        """
        if not isinstance(geom, GeometryCollection):
            return geom
        if len(geom.geoms) == 0:
            return geom
        if not all(isinstance(g, LinearRing) and g.is_valid and g.is_ring for g in geom.geoms):
            return geom
        polygons = [Polygon(ring) for ring in geom.geoms]
        return MultiPolygon(polygons)

    def _remove_z_coord_from_polygon(self,
                                     geom) -> Polygon:
        """
        Remove Z-coordinates from a Polygon.

        Args:
            geom: The Polygon geometry.

        Returns:
            Polygon: The Polygon without Z-coordinates.
        """
        geom = self._if_linear_ring_to_polygon(geom)
        exterior = [(x, y) for x, y, *_ in geom.exterior.coords]
        interiors = [[(x, y) for x, y, *_ in interior.coords]
                     for interior in geom.interiors]
        return Polygon(exterior, interiors)

    def _remove_z_coord_from_multipolygon(self,
                                          geom) -> MultiPolygon:
        """
        Remove Z-coordinates from a MultiPolygon.

        Args:
            geom: The MultiPolygon geometry.

        Returns:
            MultiPolygon: The MultiPolygon without Z-coordinates.
        """
        geom = self._if_collection_of_linear_rings_to_multipolygon(geom)
        composing_polys = [self._remove_z_coord_from_polygon(
            composing_poly) for composing_poly in geom.geoms]
        return MultiPolygon(composing_polys)

    def _remove_z_coord(self,
                        geom) -> Union[Polygon, MultiPolygon, Point, None]:
        """
        Remove Z-coordinates from a geometry.

        Args:
            geom: The geometry to process.

        Returns:
            Union[Polygon, MultiPolygon, Point, None]: The geometry without Z-coordinates.
        """
        if geom.is_empty:
            return geom
        elif isinstance(geom, (Polygon, LinearRing)):
            return self._remove_z_coord_from_polygon(geom)
        elif isinstance(geom, (MultiPolygon, GeometryCollection)):
            return self._remove_z_coord_from_multipolygon(geom)
        else:
            return geom

    def _remove_geni(self,
                     geom: Union[Polygon, MultiPolygon]) -> Union[Polygon, MultiPolygon]:
        """
        Remove "geni" (holes) from a geometry.

        Args:
            geom (Union[Polygon, MultiPolygon]): The geometry to process.

        Returns:
            Union[Polygon, MultiPolygon]: The geometry without "geni".
        """
        if geom.is_empty:
            return geom
        elif isinstance(geom, Polygon):
            return _GeniRemover(geom).trnsf_pol_all_geni()
        elif isinstance(geom, MultiPolygon):
            return MultiPolygon([_GeniRemover(pol).trnsf_pol_all_geni() for pol in geom.geoms])
        else:
            return geom

    def apply_z_coord_and_geni_removal(self,
                                       gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Apply Z-coordinate and "geni" (holes) removal to geometries in a GeoDataFrame.

        This method processes each geometry in the GeoDataFrame by removing Z-coordinates
        and holes, avoiding the deprecated apply() method that triggers warnings.

        Args:
            gdf (gpd.GeoDataFrame): The GeoDataFrame containing geometries to process.

        Returns:
            gpd.GeoDataFrame: The GeoDataFrame with processed geometries.
        """
        geometries = []
        for geom in gdf['geometry']:
            processed_geom = self._remove_z_coord(geom)
            if self.remove_geni:
                processed_geom = self._remove_geni(processed_geom)
            geometries.append(processed_geom)
        gdf['geometry'] = geometries
        return gdf


def _get_interior_point_property(self) -> gpd.GeoSeries:
    """
    Property implementation to return a point guaranteed to be inside each geometry in a GeoSeries.

    This function is designed to be used as a property getter for GeoSeries objects.
    It iterates through each geometry in the GeoSeries and generates a point that is
    guaranteed to be inside that geometry using the _get_interior_point function.

    Args:
        self: The GeoSeries instance this property is called on.

    Returns:
        gpd.GeoSeries: A new GeoSeries containing interior points for each geometry,
            preserving the original index and coordinate reference system (CRS).
    """
    result = gpd.GeoSeries(
        [_get_interior_point(geom) for geom in self],
        index=self.index,
        crs=self.crs
    )
    return result


def _get_interior_point(geometry) -> Point:
    """
    Get a point that is guaranteed to be inside a geometry.

    This function finds a point inside the given geometry, preferably the centroid
    when possible. For complex geometries where the centroid falls outside the geometry,
    it falls back to the representative_point method from Shapely.

    Algorithm:
        1. First checks if the geometry is None or empty and returns an empty Point in that case
        2. Attempts to use the centroid if it falls within the geometry (most efficient)
        3. Falls back to representative_point() if the centroid is outside the geometry
        4. In case of any exceptions, returns an empty Point

    Args:
        geometry: A Shapely geometry object (Polygon, MultiPolygon, etc.)

    Returns:
        Point: A Shapely Point that lies inside the geometry, or an empty Point if
            the geometry is None, empty, or if an error occurs.
    """
    if geometry is None or geometry.is_empty:
        return Point()
    try:
        centroid = geometry.centroid
        if geometry.contains(centroid):
            return centroid
        # If centroid is not inside, use the shapely representative_point method
        return geometry.representative_point()
    except Exception:
        try:
            return geometry.representative_point()
        except:
            return Point()


def extend_geoseries_with_interior_point() -> None:
    """
    Extends GeoPandas' GeoSeries class with an 'interior_point' property.

    This function adds a new property called 'interior_point' to the GeoPandas
    GeoSeries class if it doesn't already exist. Once added, any GeoSeries instance
    can access this property to get a series of points guaranteed to be inside
    each geometry in the series.

    The property implementation uses the _get_interior_point_property function
    as its getter, which in turn uses _get_interior_point for each geometry.

    Side Effects:
        Modifies the gpd.GeoSeries class by adding a new property if it doesn't
        already exist. This is a global change that affects all GeoSeries instances.

    Returns:
        None
    """
    if not hasattr(gpd.GeoSeries, 'interior_point'):
        gpd.GeoSeries.interior_point = property(_get_interior_point_property)


"""
Extend the GeoPandas GeoSeries class with the interior_point property.
"""
extend_geoseries_with_interior_point()


class StandardGeodataframe(gpd.GeoDataFrame):
    """
    A class to standardize and process GeoDataFrames that behaves like a GeoDataFrame itself.

    This class inherits from gpd.GeoDataFrame, allowing it to be used directly as a GeoDataFrame
    while providing additional standardization methods from _InputTransformer.

    Attributes:
        All attributes from gpd.GeoDataFrame are inherited directly
    """

    @classmethod
    def from_file(cls,
                  file_path,
                  crs="EPSG:4326",
                  remove_geni=True,
                  **kwargs):
        """
        Create a StandardGeodataframe from a file.

        Args:
            file_path (str): Path to the file to read.
            crs (str, optional): The coordinate reference system to use. Defaults to "EPSG:4326".
            remove_geni (bool, optional): Whether to remove "geni" (holes). Defaults to True.
            **kwargs: Additional arguments to pass to gpd.read_file.

        Returns:
            StandardGeodataframe: A standardized GeoDataFrame.
        """
        transformer = _InputTransformer(
            file_path, crs=crs, remove_geni=remove_geni)
        standardized_gdf = transformer.standard_gdf
        return cls(standardized_gdf)

    @classmethod
    def from_geodataframe(cls,
                          gdf,
                          crs="EPSG:4326",
                          remove_geni=True):
        """
        Create a StandardGeodataframe from an existing GeoDataFrame.

        Args:
            gdf (gpd.GeoDataFrame): The input GeoDataFrame.
            crs (str, optional): The coordinate reference system to use. Defaults to "EPSG:4326".
            remove_geni (bool, optional): Whether to remove "geni" (holes). Defaults to True.

        Returns:
            StandardGeodataframe: A standardized GeoDataFrame.
        """
        transformer = _InputTransformer(gdf, crs=crs, remove_geni=remove_geni)
        standardized_gdf = transformer.standard_gdf
        return cls(standardized_gdf)

    def __init__(self,
                 *args,
                 crs="EPSG:4326",
                 remove_geni=True,
                 **kwargs):
        """
        Initialize the StandardGeodataframe.

        This constructor either:
        1. Accepts a path or GeoDataFrame and standardizes it using _InputTransformer
        2. Passes all arguments to the parent GeoDataFrame constructor

        Args:
            *args: Arguments to pass to the parent constructor
            crs (str, optional): The coordinate reference system to use. Defaults to "EPSG:4326".
            remove_geni (bool, optional): Whether to remove "geni" (holes). Defaults to True.
            **kwargs: Keyword arguments to pass to the parent constructor
        """
        if args and (isinstance(args[0], str) or isinstance(args[0], gpd.GeoDataFrame)):
            transformer = _InputTransformer(
                args[0], crs=crs, remove_geni=remove_geni)
            super().__init__(transformer.standard_gdf, **kwargs)
        else:
            super().__init__(*args, **kwargs)

    def _validate_geometry(self,
                           geom) -> Union[Polygon, MultiPolygon]:
        """
        Validate and fix invalid geometries.

        Args:
            geom: The geometry to validate.

        Returns:
            Union[Polygon, MultiPolygon]: The validated geometry.
        """
        if not geom.is_valid:
            return geom.buffer(0)
        return geom

    def _validate_column_name(self,
                              column_name: str) -> None:
        """
        Validate the existence of a column in the GeoDataFrame.

        Args:
            column_name (str): The column name to validate.

        Raises:
            ValueError: If the column does not exist.
        """
        if column_name not in self.columns:
            raise ValueError(
                f"Column {column_name} not found in GeoDataFrame.")

    def _validate_target_value_existence(self,
                                         column_name: str,
                                         target_value: Union[str, int, float]) -> None:
        """
        Validate the existence of a target value in a column.

        Args:
            column_name (str): The column name to check.
            target_value (Union[str, int, float]): The target value to validate.

        Raises:
            ValueError: If the target value does not exist.
        """
        self._validate_column_name(column_name)
        if len(self[self[column_name] == target_value]) == 0:
            raise ValueError(
                f"No geometry found with {column_name} = {target_value}.")

    def _validate_target_value_uniqueness(self,
                                          column_name: str,
                                          target_value: Union[str, int, float]) -> None:
        """
        Validate the uniqueness of a target value in a column.

        Args:
            column_name (str): The column name to check.
            target_value (Union[str, int, float]): The target value to validate.

        Raises:
            Warning: If the target value is not unique.
        """
        self._validate_column_name(column_name)
        if len(self[self[column_name] == target_value]) > 1:
            raise Warning(
                f"More than one geometry found with {column_name} = {target_value}.")

    def _validate_target_value(self,
                               column_name: str,
                               target_value: List[Union[str, int, float]]) -> None:
        """
        Validate the existence and uniqueness of a target value in a column.

        This method ensures that the target value exists in the specified column
        and that it is unique.

        Args:
            column_name (str): The column name to check.
            target_value (List[Union[str, int, float]]): The target value to validate.

        Raises:
            ValueError: If the target value does not exist in the column.
            Warning: If the target value is not unique in the column.
        """
        self._validate_target_value_existence(column_name, target_value)
        self._validate_target_value_uniqueness(column_name, target_value)

    def _validate_geoms_sbstr_nonempty(self,
                                       column_name: str,
                                       values_to_substract: List[Union[str, int, float]],
                                       geoms_to_sbstr: Union[Polygon, MultiPolygon]) -> None:
        """
        Validate that geometries to subtract are not empty.

        This method checks if there are geometries associated with the specified
        column and values to subtract.

        Args:
            column_name (str): The column name to check.
            values_to_substract (List[Union[str, int, float]]): The values to subtract.
            geoms_to_sbstr (Union[Polygon, MultiPolygon]): The geometries to subtract.

        Raises:
            Warning: If no geometries are found for the specified column and values.
        """
        if not geoms_to_sbstr:
            raise Warning(
                f"No geometries found with {column_name} in {values_to_substract}.")

    def _get_geoms_to_substract(self,
                                column_name: str,
                                values_to_substract: List[Union[str, int, float]]) -> list:
        """
        Retrieve geometries to subtract based on column values.

        This method identifies and retrieves geometries associated with the specified
        column and values to subtract.

        Args:
            column_name (str): The column name to filter geometries.
            values_to_substract (List[Union[str, int, float]]): The values to filter geometries.

        Returns:
            list: A list containing the geometries to subtract.

        Raises:
            Warning: If no geometries are found for the specified column and values.
        """
        geometries_to_substract = []
        for value in values_to_substract:
            matches = self[self[column_name] == value]
            if len(matches) > 0:
                for idx, row in matches.iterrows():
                    geometries_to_substract.append(row.geometry)
        self._validate_geoms_sbstr_nonempty(column_name,
                                            values_to_substract,
                                            geometries_to_substract)
        return geometries_to_substract

    def _get_target_idx(self,
                        column_name: str,
                        target_value: List[Union[str, int, float]]) -> int:
        """
        Retrieve the index of the target geometry in the GeoDataFrame.

        Args:
            column_name (str): The column name to filter geometries.
            target_value (List[Union[str, int, float]]): The target value to filter geometries.

        Returns:
            int: The index of the target geometry in the GeoDataFrame.
        """
        return self[self[column_name] == target_value].index

    def _process_target_geometry(self,
                                 target_geom,
                                 column_name: str,
                                 values_to_substract: List[Union[str, int, float]]) -> Union[Polygon, MultiPolygon]:
        """
        Process the target geometry by subtracting overlapping geometries.

        This method validates the target geometry and subtracts any overlapping
        geometries associated with the specified column and values.

        Args:
            target_geom: The target geometry to process.
            column_name (str): The column name to filter geometries.
            values_to_substract (List[Union[str, int, float]]): The values to filter geometries to subtract.

        Returns:
            Union[Polygon, MultiPolygon]: The processed geometry after subtraction.
        """
        target_geom = self._validate_geometry(target_geom)
        for geom_to_sbstr in self._get_geoms_to_substract(column_name, values_to_substract):
            geom_to_sbstr = self._validate_geometry(geom_to_sbstr)
            if target_geom.intersects(geom_to_sbstr):
                target_geom = target_geom.difference(geom_to_sbstr)
        return target_geom

    def _apply_substraction(self,
                            result_gdf,
                            column_name: str,
                            target_value: List[Union[str, int, float]],
                            values_to_substract: List[Union[str, int, float]]) -> gpd.GeoDataFrame:
        """
        Apply subtraction of overlapping geometries to the target geometry.

        This method processes the target geometry by subtracting overlapping geometries
        and updates the result GeoDataFrame.

        Args:
            result_gdf (gpd.GeoDataFrame): The GeoDataFrame to update.
            column_name (str): The column name to filter geometries.
            target_value (List[Union[str, int, float]]): The target value to filter geometries.
            values_to_substract (List[Union[str, int, float]]): The values to filter geometries to subtract.

        Returns:
            gpd.GeoDataFrame: The updated GeoDataFrame after subtraction.
        """
        target_indices = self._get_target_idx(column_name, target_value)
        for idx in target_indices:
            result_gdf.loc[idx, 'geometry'] = self._process_target_geometry(
                result_gdf.loc[idx, 'geometry'],
                column_name,
                values_to_substract
            )
        return result_gdf

    def _calculate_surface_area(self,
                                conversion_factor=1.0,
                                column_name='SURF_A') -> Tuple['StandardGeodataframe', str]:
        """
        Helper method to calculate surface area with conversion factor.

        This method projects the geometries to EPSG:3395 (World Mercator projection)
        to calculate accurate areas, then applies the conversion factor and
        returns the geometries to their original projection.

        Args:
            conversion_factor (float): Factor to convert square meters to desired units.
                Use 1.0 for square meters, 10^6 for square kilometers.
            column_name (str): Name of the output column to store the calculated areas.

        Returns:
            tuple: (StandardGeodataframe, str) - A tuple containing:
                - The processed GeoDataFrame with the new area column
                - The original CRS of the GeoDataFrame
        """
        original_crs = self.crs
        result = self.copy()
        result = result.to_crs(epsg=3395)
        result[column_name] = result['geometry'].area / conversion_factor
        result = result.to_crs(original_crs)
        return result, original_crs

    def get_interior_points(self,
                            geometry_column='geometry',
                            inplace=False,
                            column_name='interior_point') -> Union[gpd.GeoSeries, 'StandardGeodataframe']:
        """
        Calculate interior points for each geometry.

        Args:
            geometry_column (str): Name of the geometry column to process. Defaults to 'geometry'.
            inplace (bool): If True, adds the interior points as a new column to the GeoDataFrame.
                If False, returns the interior points as a GeoSeries. Defaults to False.
            column_name (str): Name of the column to store interior points when inplace=True.
                Defaults to 'interior_point'.

        Returns:
            Union[gpd.GeoSeries, 'StandardGeodataframe']: 
                If inplace=False, returns a GeoSeries containing interior points.
                If inplace=True, returns self with a new column containing interior points.

        Raises:
            ValueError: If the geometry_column doesn't exist in the GeoDataFrame.
        """
        self._validate_column_name(geometry_column)
        interior_points = self[geometry_column].interior_point
        if inplace:
            self[column_name] = interior_points
            return self
        else:
            return interior_points

    def substract_overlapping_geometries(self,
                                         column_name: str,
                                         args: Union[Tuple[List[Union[str, int, float]],
                                                     List[Union[str, int, float]]],
                                                     dict],
                                         remove_geni: bool = True,
                                         inplace: bool = False) -> 'StandardGeodataframe':
        """
        Subtract overlapping geometries from target geometries.

        This method subtracts specified geometries (holes or overlapping areas) from 
        target geometries in the GeoDataFrame. The subtraction can be specified using 
        either a tuple or a dictionary.

        Args:
            column_name (str): The column name to identify geometries.
            args (Union[Tuple[List[Union[str, int, float]], List[Union[str, int, float]]], dict]): 
                A tuple or dictionary specifying target and subtracting geometries.
                - Tuple: A pair of lists where the first list contains target values 
                and the second list contains values to subtract.
                - Dictionary: A mapping where keys are target values and values are 
                lists of geometries to subtract.
            remove_geni (bool, optional): Whether to remove holes from resulting geometries.
                Defaults to True.
            inplace (bool, optional): If True, modifies the current GeoDataFrame.
                If False, returns a new GeoDataFrame. Defaults to False.

        Returns:
            StandardGeodataframe: If inplace=False, returns a new instance with subtracted geometries.
                                If inplace=True, returns self with updated geometries.

        Raises:
            ValueError: If the `args` format is invalid or if column_name doesn't exist.
        """
        self._validate_column_name(column_name)
        result_gdf = self.copy()
        if isinstance(args, dict):
            errors = []
            for target_value, values_to_substract in args.items():
                try:
                    result_gdf = self._apply_substraction(
                        result_gdf, column_name, target_value, values_to_substract)
                except Exception as e:
                    errors.append(
                        f"Error processing {target_value} with {values_to_substract}: {str(e)}")
            if errors:
                print(f"Encountered {len(errors)} errors during processing:")
                for error in errors:
                    print(f"  - {error}")
        elif isinstance(args, tuple) and len(args) == 2:
            target_value, values_to_substract = args
            result_gdf = self._apply_substraction(
                result_gdf, column_name, target_value, values_to_substract)
        else:
            raise ValueError(
                "Invalid args format. Must be a tuple (target_value, values_to_substract) or a dictionary.")
        if inplace:
            self.geometry = result_gdf.geometry
            return self
        else:
            return StandardGeodataframe(result_gdf, crs=self.crs, remove_geni=remove_geni)

    def get_m2_surface_area(self,
                            inplace=False,
                            column_name='SURF_A_M2'):
        """
        Calculate surface area in square meters for each geometry.

        This method calculates the surface area of each geometry in the GeoDataFrame in square meters 
        using the EPSG:3395 CRS (World Mercator projection).

        Args:
            inplace (bool): If True, adds the surface area as a new column to the GeoDataFrame.
                If False, returns a copy of the GeoDataFrame with the added column. Defaults to False.
            column_name (str): Name of the column to store the surface area. Defaults to 'SURF_A_M2'.

        Returns:
            StandardGeodataframe: If inplace=False, a copy of the GeoDataFrame with the added column.
                                If inplace=True, self with the added column.
        """
        result, _ = self._calculate_surface_area(1.0, column_name)
        if inplace:
            self[column_name] = result[column_name]
            return self
        else:
            return StandardGeodataframe(result, crs=self.crs)

    def get_km2_surface_area(self,
                             inplace=False,
                             column_name='SURF_A_KM2'):
        """
        Calculate surface area in square kilometers for each geometry.

        This method calculates the surface area of each geometry in the GeoDataFrame in square meters 
        using the EPSG:3395 CRS (World Mercator projection) and converts it to square kilometers.

        Args:
            inplace (bool): If True, adds the surface area as a new column to the GeoDataFrame.
                If False, returns a copy of the GeoDataFrame with the added column. Defaults to False.
            column_name (str): Name of the column to store the surface area. Defaults to 'SURF_A_KM2'.

        Returns:
            StandardGeodataframe: If inplace=False, a copy of the GeoDataFrame with the added column.
                                If inplace=True, self with the added column.
        """
        result, _ = self._calculate_surface_area(10**6, column_name)
        if inplace:
            self[column_name] = result[column_name]
            return self
        else:
            return StandardGeodataframe(result, crs=self.crs)
