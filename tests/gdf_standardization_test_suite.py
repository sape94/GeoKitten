"""
Test suite for the gdf_standardization module.

This module contains comprehensive unit tests for all classes, methods, and functions
in the gdf_standardization module. The tests cover functionality related to
GeoDataFrame standardization, geometry processing (Z-coordinate and hole removal),
calculating surface areas, finding interior points, and subtracting overlapping geometries.

Test Classes:
    TestGlobalFunctions: Tests for standalone functions in gdf_standardization
    TestGeoDataFrameAdapter: Tests for the GeoDataFrameAdapter class
    TestGeniRemover: Tests for the _GeniRemover class
    TestInputTransformer: Tests for the _InputTransformer class
    TestStandardGeodataframe: Tests for the StandardGeodataframe class

Each test class includes setup methods and multiple test cases that validate
the functionality of the corresponding module components. Tests use sample geometries
and test files to verify correct behavior for standard use cases and edge cases.

Author: Sergio A. Pelayo Escalera
Created: 2025-04-09
Last-modified: 2025-04-09
"""

from geokitten.gdf_standardization import (
    GeoDataFrameAdapter,
    _GeniRemover,
    _InputTransformer,
    _get_interior_point,
    extend_geoseries_with_interior_point,
    StandardGeodataframe
)
import unittest
import os
import sys
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon, LinearRing, GeometryCollection
from pathlib import Path

# Add the src directory to the Python path
sys.path.append(str(Path(__file__).parent.parent))


class TestGlobalFunctions(unittest.TestCase):
    """
    Test suite for global functions in gdf_standardization.py.

    This test class verifies the functionality of standalone functions in the module,
    including _get_interior_point and extend_geoseries_with_interior_point.

    Attributes:
        polygon (Polygon): A simple square polygon for testing.
        empty_geometry (Polygon): An empty polygon for edge case testing.
    """

    def setUp(self):
        """
        Set up test fixtures before each test.

        Creates geometric objects used across multiple test methods:
        - A simple square polygon
        - An empty polygon
        """
        self.polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
        self.empty_geometry = Polygon()

    def test_get_interior_point(self):
        """
        Test _get_interior_point function.

        Verifies that:
        1. For a valid polygon, the function returns a Point at its centroid
        2. For an empty geometry, the function returns an empty Point
        """
        # Test with valid polygon - should return centroid
        point = _get_interior_point(self.polygon)
        self.assertIsInstance(point, Point)
        self.assertEqual(point.x, 0.5)
        self.assertEqual(point.y, 0.5)

        # Test with empty geometry - should return empty Point
        point = _get_interior_point(self.empty_geometry)
        self.assertIsInstance(point, Point)
        self.assertTrue(point.is_empty)

    def test_extend_geoseries_with_interior_point(self):
        """
        Test that extend_geoseries_with_interior_point adds the property.

        Verifies that:
        1. The function adds an 'interior_point' property to the gpd.GeoSeries class
        2. The property returns a GeoSeries of interior points when accessed
        3. The returned points are valid Shapely Point objects
        """
        extend_geoseries_with_interior_point()
        self.assertTrue(hasattr(gpd.GeoSeries, 'interior_point'))

        # Create a GeoSeries and test the property
        gs = gpd.GeoSeries([self.polygon])
        self.assertIsInstance(gs.interior_point, gpd.GeoSeries)
        self.assertEqual(len(gs.interior_point), 1)
        self.assertIsInstance(gs.interior_point.iloc[0], Point)


class TestGeoDataFrameAdapter(unittest.TestCase):
    """
    Test suite for GeoDataFrameAdapter class.

    This test class verifies the functionality of the GeoDataFrameAdapter class,
    which provides a unified interface for working with GeoPandas GeoDataFrames.

    Attributes:
        test_file_path (str): Path to the test shapefile.
        gdf (gpd.GeoDataFrame): A test GeoDataFrame loaded from the shapefile.
        adapter (GeoDataFrameAdapter): An adapter instance wrapping the test GeoDataFrame.
    """

    def setUp(self):
        """
        Set up test fixtures before each test.

        Initializes:
        1. Path to a test shapefile
        2. A GeoDataFrame loaded from the test file
        3. A GeoDataFrameAdapter instance wrapping the GeoDataFrame
        """
        current_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_file_path = os.path.join(
            current_directory,
            'tests_files',
            'inputs',
            'gdf_standardization_test_file.shp'
        )
        self.gdf = gpd.read_file(self.test_file_path)
        self.adapter = GeoDataFrameAdapter(self.gdf)

    def test_initialization(self):
        """
        Test initialization of GeoDataFrameAdapter.

        Verifies that:
        1. The framework property is set correctly
        2. The adapter properly wraps the GeoDataFrame
        3. The wrapped GeoDataFrame has identical properties to the original
        """
        self.assertEqual(self.adapter.framework, "geopandas")

        # Instead of comparing the entire DataFrame directly, check specific attributes
        self.assertTrue(self.adapter.geodataframe.equals(self.gdf))
        self.assertEqual(len(self.adapter.geodataframe), len(self.gdf))
        self.assertEqual(list(self.adapter.geodataframe.columns),
                         list(self.gdf.columns))

        # Check that the CRS matches
        self.assertEqual(self.adapter.geodataframe.crs, self.gdf.crs)

        # Verify object identity (they are the same object)
        self.assertIs(self.adapter.geodataframe, self.gdf)

    def test_initialization_with_invalid_input(self):
        """
        Test initialization with invalid input.

        Verifies that attempting to initialize with a non-GeoDataFrame
        object raises a TypeError.
        """
        with self.assertRaises(TypeError):
            GeoDataFrameAdapter("not_a_geodataframe")

    def test_to_crs(self):
        """
        Test to_crs method.

        Verifies that:
        1. Without inplace=True: Returns a transformed GeoDataFrame with the new CRS
        2. With inplace=True: Modifies the original GeoDataFrame and returns it
        """
        # Test without inplace
        result = self.adapter.to_crs("EPSG:3857")
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:3857")

        # Test with inplace
        original_crs = self.adapter.crs
        result = self.adapter.to_crs("EPSG:4326", inplace=True)
        self.assertEqual(self.adapter.crs, "EPSG:4326")
        self.assertTrue(result.equals(self.adapter.geodataframe))

    def test_crs_property(self):
        """
        Test crs property.

        Verifies that:
        1. The getter returns the correct CRS of the GeoDataFrame
        2. The setter correctly changes the CRS of the GeoDataFrame
        """
        self.assertEqual(self.adapter.crs, self.gdf.crs)

        # Test setter
        self.adapter.crs = "EPSG:3857"
        self.assertEqual(self.adapter.crs, "EPSG:3857")

    def test_copy(self):
        """
        Test copy method.

        Verifies that the method returns a new GeoDataFrame that is a deep copy
        of the original, with identical data but different object identity.
        """
        copy = self.adapter.copy()
        self.assertIsInstance(copy, gpd.GeoDataFrame)
        self.assertIsNot(copy, self.adapter.geodataframe)

    def test_columns_property(self):
        """
        Test columns property.

        Verifies that the property returns a list of the column names
        in the GeoDataFrame.
        """
        self.assertEqual(self.adapter.columns, list(self.gdf.columns))

    def test_getitem(self):
        """
        Test __getitem__ method.

        Verifies that:
        1. With a column name: Returns the corresponding column as a Series/GeoSeries
        2. With a slice: Returns a subset of the GeoDataFrame
        """
        # Test with column name
        result = self.adapter["geometry"]
        self.assertIsInstance(result, gpd.GeoSeries)

        # Test with slice
        result = self.adapter[0:2]
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(len(result), 2)

    def test_iterrows(self):
        """
        Test iterrows method.

        Verifies that the method correctly iterates over rows of the GeoDataFrame,
        yielding (index, Series) pairs for each row.
        """
        count = 0
        for idx, row in self.adapter.iterrows():
            self.assertIsInstance(row, pd.Series)
            count += 1
        self.assertEqual(count, len(self.gdf))

    def test_apply(self):
        """
        Test apply method.

        Verifies that:
        1. The method applies a function along an axis of the GeoDataFrame
        2. Works correctly with both axis=0 (columns) and axis=1 (rows)
        """
        # Simple function to apply
        def test_func(x):
            return len(x)

        # Apply to columns
        result = self.adapter.apply(test_func)
        self.assertIsInstance(result, pd.Series)

        # Apply to rows
        result = self.adapter.apply(test_func, axis=1)
        self.assertIsInstance(result, pd.Series)

    def test_to_native(self):
        """
        Test to_native method.

        Verifies that the method returns the underlying GeoDataFrame object.
        """
        native = self.adapter.to_native()
        self.assertIs(native, self.gdf)

    def test_from_file(self):
        """
        Test from_file class method.

        Verifies that the method creates a GeoDataFrameAdapter instance
        by reading a GeoDataFrame from a file path.
        """
        adapter = GeoDataFrameAdapter.from_file(self.test_file_path)
        self.assertIsInstance(adapter, GeoDataFrameAdapter)
        self.assertEqual(adapter.framework, "geopandas")

    def test_from_file_with_invalid_path(self):
        """
        Test from_file with invalid path.

        Verifies that attempting to create an adapter from a non-existent
        file path raises a ValueError.
        """
        with self.assertRaises(ValueError):
            GeoDataFrameAdapter.from_file("not_a_file.shp")


class TestGeniRemover(unittest.TestCase):
    """
    Test suite for _GeniRemover class.

    This test class verifies the functionality of the _GeniRemover class,
    which is responsible for removing holes from polygons.

    Attributes:
        polygon_with_hole (Polygon): A test polygon containing a hole.
        remover (_GeniRemover): A _GeniRemover instance for the test polygon.
        simple_polygon (Polygon): A simple polygon without holes.
    """

    def setUp(self):
        """
        Set up test fixtures before each test.

        Creates:
        1. A polygon with a hole
        2. A _GeniRemover instance for the polygon with a hole
        3. A simple polygon without holes
        """
        # Create a polygon with a hole
        exterior = [(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)]
        hole = [(3, 3), (3, 7), (7, 7), (7, 3), (3, 3)]
        self.polygon_with_hole = Polygon(exterior, [hole])
        self.remover = _GeniRemover(self.polygon_with_hole)

        # Simple polygon without holes
        self.simple_polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])

    def test_initialization(self):
        """
        Test initialization of _GeniRemover.

        Verifies that the geometry is correctly stored in the instance.
        """
        self.assertEqual(self.remover.geom, self.polygon_with_hole)

    def test_calculate_distance(self):
        """
        Test _calculate_distance method.

        Verifies that the method correctly calculates the Euclidean distance
        between two points.
        """
        distance = self.remover._calculate_distance((0, 0), (3, 4))
        self.assertEqual(distance, 5.0)

    def test_initialize_nearest_search(self):
        """
        Test _initialize_nearest_search method.

        Verifies that the method initializes the nearest pair search with
        null results and infinite distance.
        """
        result, dist = self.remover._initialize_nearest_search()
        self.assertIsNone(result)
        self.assertEqual(dist, float('inf'))

    def test_update_nearest_pair(self):
        """
        Test _update_nearest_pair method.

        Verifies that:
        1. The method updates the nearest pair when a closer pair is found
        2. The method does not update when the new pair is farther
        """
        nearest_pair, min_dist = ((0, 0), (3, 3)), 5.0

        # This pair is closer, should update
        new_pair, new_dist = self.remover._update_nearest_pair(
            (1, 1), (3, 3), nearest_pair, min_dist
        )
        self.assertEqual(new_pair, ((1, 1), (3, 3)))
        self.assertLess(new_dist, min_dist)

        # This pair is farther, should not update
        new_pair, new_dist = self.remover._update_nearest_pair(
            (0, 0), (10, 10), nearest_pair, min_dist
        )
        self.assertEqual(new_pair, nearest_pair)
        self.assertEqual(new_dist, min_dist)

    def test_find_nearest_points(self):
        """
        Test _find_nearest_points method.

        Verifies that the method finds the closest pair of points between
        a polygon's exterior and a hole.
        """
        exterior = list(self.polygon_with_hole.exterior.coords)
        hole = list(self.polygon_with_hole.interiors[0].coords)

        ext_point, hole_point, dist = self.remover._find_nearest_points(
            exterior, hole)

        self.assertIsInstance(ext_point, tuple)
        self.assertIsInstance(hole_point, tuple)
        self.assertIsInstance(dist, float)
        self.assertGreaterEqual(dist, 0.0)

    def test_trnsf_pol_all_geni(self):
        """
        Test trnsf_pol_all_geni method.

        Verifies that:
        1. For a polygon with a hole, the method returns a polygon without holes
        2. For a simple polygon, the method returns the same polygon unchanged
        """
        # With a polygon that has a hole
        result = self.remover.trnsf_pol_all_geni()
        self.assertIsInstance(result, Polygon)
        self.assertEqual(len(result.interiors), 0)  # No more holes

        # With a simple polygon without holes
        remover = _GeniRemover(self.simple_polygon)
        result = remover.trnsf_pol_all_geni()
        self.assertIsInstance(result, Polygon)
        self.assertEqual(len(result.interiors), 0)
        self.assertEqual(result, self.simple_polygon)  # Should be unchanged


class TestInputTransformer(unittest.TestCase):
    """
    Test suite for _InputTransformer class.

    This test class verifies the functionality of the _InputTransformer class,
    which transforms input data (file paths or GeoDataFrames) into standardized
    GeoDataFrames.

    Attributes:
        test_file_path (str): Path to the test shapefile.
        gdf (gpd.GeoDataFrame): A test GeoDataFrame loaded from the shapefile.
    """

    def setUp(self):
        """
        Set up test fixtures before each test.

        Initializes:
        1. Path to a test shapefile
        2. A GeoDataFrame loaded from the test file
        """
        current_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_file_path = os.path.join(
            current_directory,
            'tests_files',
            'inputs',
            'gdf_standardization_test_file.shp'
        )
        self.gdf = gpd.read_file(self.test_file_path)

    def test_initialization_from_file(self):
        """
        Test initialization from a file path.

        Verifies that the transformer correctly creates a GeoDataFrame
        from a file path.
        """
        transformer = _InputTransformer(self.test_file_path)
        self.assertIsInstance(transformer.standard_gdf, gpd.GeoDataFrame)

    def test_initialization_from_geodataframe(self):
        """
        Test initialization from a GeoDataFrame.

        Verifies that the transformer correctly accepts and processes
        an existing GeoDataFrame.
        """
        transformer = _InputTransformer(self.gdf)
        self.assertIsInstance(transformer.standard_gdf, gpd.GeoDataFrame)

    def test_initialization_with_invalid_input(self):
        """
        Test initialization with invalid input.

        Verifies that attempting to initialize with an invalid input type
        raises a ValueError.
        """
        with self.assertRaises(ValueError):
            _InputTransformer(123)

    def test_is_path(self):
        """
        Test _is_path method.

        Verifies that the method correctly identifies whether the input
        is a file path.
        """
        transformer = _InputTransformer(self.test_file_path)
        self.assertTrue(transformer._is_path())

        transformer = _InputTransformer(self.gdf)
        self.assertFalse(transformer._is_path())

    def test_is_geodataframe(self):
        """
        Test _is_geodataframe method.

        Verifies that the method correctly identifies whether the input
        is a GeoDataFrame.
        """
        transformer = _InputTransformer(self.gdf)
        self.assertTrue(transformer._is_geodataframe())

        transformer = _InputTransformer(self.test_file_path)
        self.assertFalse(transformer._is_geodataframe())

    def test_path_validation(self):
        """
        Test _path_validation method.

        Verifies that:
        1. The method does not raise exceptions for valid file paths
        2. The method raises ValueError for non-existent paths
        """
        transformer = _InputTransformer(self.test_file_path)
        # Should not raise an exception
        transformer._path_validation()

        # Should raise ValueError for non-existent path
        transformer.gdf_input = "non_existent_file.shp"
        with self.assertRaises(ValueError):
            transformer._path_validation()

    def test_get_gdf(self):
        """
        Test _get_gdf method.

        Verifies that the method returns a GeoDataFrame regardless of whether
        the input was a file path or a GeoDataFrame.
        """
        # From file path
        transformer = _InputTransformer(self.test_file_path)
        gdf = transformer._get_gdf()
        self.assertIsInstance(gdf, gpd.GeoDataFrame)

        # From GeoDataFrame
        transformer = _InputTransformer(self.gdf)
        gdf = transformer._get_gdf()
        self.assertIsInstance(gdf, gpd.GeoDataFrame)
        self.assertIs(gdf, self.gdf)

    def test_set_standard_crs(self):
        """
        Test set_standard_crs method.

        Verifies that:
        1. With default CRS: Sets the CRS to EPSG:4326
        2. With custom CRS: Sets the CRS to the specified value
        """
        # With default CRS
        transformer = _InputTransformer(self.gdf)
        gdf = transformer.set_standard_crs()
        self.assertEqual(gdf.crs, "EPSG:4326")

        # With custom CRS
        transformer = _InputTransformer(self.gdf, crs="EPSG:3857")
        gdf = transformer.set_standard_crs()
        self.assertEqual(gdf.crs, "EPSG:3857")

    def test_if_linear_ring_to_polygon(self):
        """
        Test _if_linear_ring_to_polygon method.

        Verifies that:
        1. When given a LinearRing, the method converts it to a Polygon
        2. When given other geometry types, the method returns them unchanged
        """
        transformer = _InputTransformer(self.gdf)

        # Test with LinearRing
        ring = LinearRing([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
        result = transformer._if_linear_ring_to_polygon(ring)
        self.assertIsInstance(result, Polygon)

        # Test with non-LinearRing
        point = Point(0, 0)
        result = transformer._if_linear_ring_to_polygon(point)
        self.assertEqual(result, point)

    def test_if_collection_of_linear_rings_to_multipolygon(self):
        """
        Test _if_collection_of_linear_rings_to_multipolygon method.

        Verifies that:
        1. When given a GeometryCollection of LinearRings, converts to MultiPolygon
        2. When given other geometry types, returns them unchanged
        """
        transformer = _InputTransformer(self.gdf)

        # Test with GeometryCollection of LinearRings
        rings = [
            LinearRing([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]),
            LinearRing([(2, 2), (2, 3), (3, 3), (3, 2), (2, 2)])
        ]
        collection = GeometryCollection(rings)
        result = transformer._if_collection_of_linear_rings_to_multipolygon(
            collection)
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)

        # Test with non-GeometryCollection
        point = Point(0, 0)
        result = transformer._if_collection_of_linear_rings_to_multipolygon(
            point)
        self.assertEqual(result, point)

    def test_remove_z_coord_from_polygon(self):
        """
        Test _remove_z_coord_from_polygon method.

        Verifies that the method removes Z coordinates from a Polygon,
        returning a new Polygon with only X and Y coordinates.
        """
        transformer = _InputTransformer(self.gdf)

        # Polygon with Z coordinates
        polygon_with_z = Polygon(
            [(0, 0, 1), (0, 1, 1), (1, 1, 1), (1, 0, 1), (0, 0, 1)])
        result = transformer._remove_z_coord_from_polygon(polygon_with_z)

        # Check result is a Polygon without Z coordinates
        self.assertIsInstance(result, Polygon)
        for coord in result.exterior.coords:
            self.assertEqual(len(coord), 2)

    def test_remove_z_coord_from_multipolygon(self):
        """
        Test _remove_z_coord_from_multipolygon method.

        Verifies that the method removes Z coordinates from all component
        Polygons in a MultiPolygon.
        """
        transformer = _InputTransformer(self.gdf)

        # MultiPolygon with Z coordinates
        polygons_with_z = [
            Polygon([(0, 0, 1), (0, 1, 1), (1, 1, 1), (1, 0, 1), (0, 0, 1)]),
            Polygon([(2, 2, 1), (2, 3, 1), (3, 3, 1), (3, 2, 1), (2, 2, 1)])
        ]
        multi_polygon = MultiPolygon(polygons_with_z)
        result = transformer._remove_z_coord_from_multipolygon(multi_polygon)

        # Check result is a MultiPolygon without Z coordinates
        self.assertIsInstance(result, MultiPolygon)
        for polygon in result.geoms:
            for coord in polygon.exterior.coords:
                self.assertEqual(len(coord), 2)

    def test_remove_z_coord(self):
        """
        Test _remove_z_coord method.

        Verifies that:
        1. For empty geometries: Returns an empty geometry
        2. For Points: Returns the original Point (not processed)
        3. For Polygons: Removes Z coordinates
        4. For MultiPolygons: Removes Z coordinates from all Polygons
        """
        transformer = _InputTransformer(self.gdf)

        # Test with empty geometry
        empty = Polygon()
        result = transformer._remove_z_coord(empty)
        self.assertTrue(result.is_empty)

        # Test with Point
        point = Point(0, 0, 1)
        result = transformer._remove_z_coord(point)
        self.assertEqual(result, point)  # Points are not processed

        # Test with Polygon
        polygon = Polygon(
            [(0, 0, 1), (0, 1, 1), (1, 1, 1), (1, 0, 1), (0, 0, 1)])
        result = transformer._remove_z_coord(polygon)
        self.assertIsInstance(result, Polygon)

        # Test with MultiPolygon
        multi_polygon = MultiPolygon([polygon])
        result = transformer._remove_z_coord(multi_polygon)
        self.assertIsInstance(result, MultiPolygon)

    def test_remove_geni(self):
        """
        Test _remove_geni method.

        Verifies that:
        1. For empty geometries: Returns an empty geometry
        2. For simple Polygons: Returns the original Polygon unchanged
        3. For Polygons with holes: Removes all holes
        4. For MultiPolygons: Removes holes from all component Polygons
        """
        transformer = _InputTransformer(self.gdf)

        # Test with empty geometry
        empty = Polygon()
        result = transformer._remove_geni(empty)
        self.assertTrue(result.is_empty)

        # Test with simple Polygon
        simple = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
        result = transformer._remove_geni(simple)
        self.assertEqual(result, simple)  # Should be unchanged

        # Test with Polygon with hole
        exterior = [(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)]
        hole = [(3, 3), (3, 7), (7, 7), (7, 3), (3, 3)]
        with_hole = Polygon(exterior, [hole])
        result = transformer._remove_geni(with_hole)
        self.assertEqual(len(result.interiors), 0)  # No more holes

        # Test with MultiPolygon
        multi = MultiPolygon([with_hole, simple])
        result = transformer._remove_geni(multi)
        self.assertIsInstance(result, MultiPolygon)
        for polygon in result.geoms:
            self.assertEqual(len(polygon.interiors), 0)  # No holes

    def test_apply_z_coord_and_geni_removal_with_geni_removal(self):
        """
        Test apply_z_coord_and_geni_removal method with remove_geni=True.

        Verifies that when remove_geni=True, the method:
        1. Removes Z coordinates from all geometries
        2. Removes holes from all Polygons
        """
        transformer = _InputTransformer(self.gdf)

        # Create test GeoDataFrame with Z coords and holes
        exterior = [(0, 0, 1), (0, 10, 1), (10, 10, 1), (10, 0, 1), (0, 0, 1)]
        hole = [(3, 3, 1), (3, 7, 1), (7, 7, 1), (7, 3, 1), (3, 3, 1)]
        with_hole = Polygon(exterior, [hole])

        simple = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])

        test_gdf = gpd.GeoDataFrame(geometry=[with_hole, simple])

        # Test with remove_geni=True
        transformer.remove_geni = True
        result = transformer.apply_z_coord_and_geni_removal(test_gdf)

        self.assertIsInstance(result, gpd.GeoDataFrame)
        for geom in result['geometry']:
            if isinstance(geom, Polygon):
                self.assertEqual(len(geom.interiors), 0)  # No holes
            for coord in geom.exterior.coords:
                self.assertEqual(len(coord), 2)  # No Z coords

    def test_apply_z_coord_and_geni_removal_without_geni_removal(self):
        """
        Test apply_z_coord_and_geni_removal method with remove_geni=False.

        Verifies that when remove_geni=False, the method:
        1. Removes Z coordinates from all geometries
        2. Preserves holes in Polygons
        """
        transformer = _InputTransformer(self.gdf)

        # Create test GeoDataFrame with Z coords and holes
        exterior = [(0, 0, 1), (0, 10, 1), (10, 10, 1), (10, 0, 1), (0, 0, 1)]
        hole = [(3, 3, 1), (3, 7, 1), (7, 7, 1), (7, 3, 1), (3, 3, 1)]
        with_hole = Polygon(exterior, [hole])

        simple = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])

        test_gdf = gpd.GeoDataFrame(geometry=[with_hole, simple])

        # Test with remove_geni=False
        transformer.remove_geni = False
        result = transformer.apply_z_coord_and_geni_removal(test_gdf)

        self.assertIsInstance(result, gpd.GeoDataFrame)
        # The first geometry should still have a hole
        self.assertEqual(len(result['geometry'].iloc[0].interiors), 1)
        # But no Z coordinates
        for geom in result['geometry']:
            for coord in geom.exterior.coords:
                self.assertEqual(len(coord), 2)


class TestStandardGeodataframe(unittest.TestCase):
    """
    Test suite for StandardGeodataframe class.

    This test class verifies the functionality of the StandardGeodataframe class,
    which extends GeoPandas GeoDataFrame with additional functionality for
    standardized geographic data processing.

    Attributes:
        test_file_path (str): Path to the test shapefile.
        gdf (StandardGeodataframe): A test StandardGeodataframe loaded from the shapefile.
        simple_polygon (Polygon): A simple polygon for testing geometry validation.
    """

    def setUp(self):
        """
        Set up test fixtures before each test method.

        Initializes:
        1. Path to a test shapefile
        2. A StandardGeodataframe loaded from the test file
        3. A simple polygon for testing
        """
        # Crear path del archivo de prueba
        current_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_file_path = os.path.join(
            current_directory,
            'tests_files',
            'inputs',
            'gdf_standardization_test_file.shp'
        )

        # Cargar datos para las pruebas
        self.gdf = StandardGeodataframe(self.test_file_path)

        # Crear un polígono simple para pruebas
        self.simple_polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])

    def test_initialization_from_file(self):
        """
        Test initialization from a shapefile.

        Verifies that:
        1. The object is an instance of StandardGeodataframe
        2. The CRS is set to the default EPSG:4326
        """
        self.assertIsInstance(self.gdf, StandardGeodataframe)
        self.assertEqual(self.gdf.crs, "EPSG:4326")

    def test_initialization_from_geodataframe(self):
        """
        Test initialization from an existing GeoDataFrame.

        Verifies that a StandardGeodataframe can be created from
        an existing GeoPandas GeoDataFrame.
        """
        basic_gdf = gpd.read_file(self.test_file_path)
        standard_gdf = StandardGeodataframe(basic_gdf)
        self.assertIsInstance(standard_gdf, StandardGeodataframe)

    def test_from_file_class_method(self):
        """
        Test from_file class method.

        Verifies that the method creates a StandardGeodataframe from
        a file path.
        """
        gdf = StandardGeodataframe.from_file(self.test_file_path)
        self.assertIsInstance(gdf, StandardGeodataframe)
        self.assertEqual(gdf.crs, "EPSG:4326")

    def test_from_geodataframe_class_method(self):
        """
        Test from_geodataframe class method.

        Verifies that the method creates a StandardGeodataframe from
        an existing GeoDataFrame.
        """
        basic_gdf = gpd.read_file(self.test_file_path)
        gdf = StandardGeodataframe.from_geodataframe(basic_gdf)
        self.assertIsInstance(gdf, StandardGeodataframe)

    def test_validate_geometry(self):
        """
        Test _validate_geometry method.

        Verifies that:
        1. Valid geometries are returned unchanged
        2. Invalid geometries (e.g., self-intersecting) are corrected
        """
        # Valid geometry should be returned unchanged
        result = self.gdf._validate_geometry(self.simple_polygon)
        self.assertEqual(result, self.simple_polygon)

        # Create an invalid geometry (self-intersecting)
        invalid = Polygon([(0, 0), (1, 1), (0, 1), (1, 0)])
        self.assertFalse(invalid.is_valid)

        # _validate_geometry should fix it
        result = self.gdf._validate_geometry(invalid)
        self.assertTrue(result.is_valid)

    def test_validate_column_name(self):
        """
        Test _validate_column_name method.

        Verifies that:
        1. Valid column names do not raise exceptions
        2. Non-existent column names raise ValueError
        """
        # Valid column should not raise exception
        existing_column = self.gdf.columns[0]
        try:
            self.gdf._validate_column_name(existing_column)
            validation_passed = True
        except ValueError:
            validation_passed = False
        self.assertTrue(validation_passed)

        # Invalid column should raise ValueError
        with self.assertRaises(ValueError):
            self.gdf._validate_column_name("non_existent_column")

    def test_validate_target_value_existence(self):
        """
        Test _validate_target_value_existence method.

        Verifies that:
        1. Existing values in a column do not raise exceptions
        2. Non-existent values raise ValueError
        """
        # Use a column that exists
        column_name = "Name"  # Assuming this exists in the test data

        # Get an existing value
        existing_value = self.gdf[column_name].iloc[0]

        # Valid value should not raise exception
        try:
            self.gdf._validate_target_value_existence(
                column_name, existing_value)
            validation_passed = True
        except ValueError:
            validation_passed = False
        self.assertTrue(validation_passed)

        # Invalid value should raise ValueError
        with self.assertRaises(ValueError):
            self.gdf._validate_target_value_existence(
                column_name, "non_existent_value")

    def test_get_interior_points_without_inplace(self):
        """
        Test get_interior_points method without inplace option.

        Verifies that with inplace=False, the method returns a GeoSeries
        containing interior points for each geometry.
        """
        interior_points = self.gdf.get_interior_points()
        self.assertIsInstance(interior_points, gpd.GeoSeries)
        self.assertEqual(len(interior_points), len(self.gdf))

    def test_get_interior_points_with_inplace(self):
        """
        Test get_interior_points method with inplace=True.

        Verifies that with inplace=True, the method adds a new column to
        the GeoDataFrame containing interior points.
        """
        result = self.gdf.get_interior_points(
            inplace=True, column_name='test_points')
        self.assertIsInstance(result, StandardGeodataframe)
        self.assertIn('test_points', result.columns)

    def test_get_interior_points_with_invalid_column(self):
        """
        Test get_interior_points with invalid column.

        Verifies that using a non-existent column name raises ValueError.
        """
        with self.assertRaises(ValueError):
            self.gdf.get_interior_points(geometry_column="non_existent_column")

    def test_get_m2_surface_area_without_inplace(self):
        """
        Test get_m2_surface_area method without inplace option.

        Verifies that with inplace=False, the method returns a new GeoDataFrame
        with an additional column containing area values in square meters.
        """
        result = self.gdf.get_m2_surface_area()
        self.assertIsInstance(result, StandardGeodataframe)
        self.assertIn('SURF_A_M2', result.columns)

        # Verify values are positive
        self.assertTrue((result['SURF_A_M2'] > 0).all())

    def test_get_m2_surface_area_with_inplace(self):
        """
        Test get_m2_surface_area method with inplace=True.

        Verifies that with inplace=True, the method adds a new column to
        the GeoDataFrame containing area values in square meters.
        """
        original_columns = self.gdf.columns.tolist()
        result = self.gdf.get_m2_surface_area(
            inplace=True, column_name='area_test')
        self.assertIsInstance(result, StandardGeodataframe)
        self.assertIn('area_test', result.columns)
        self.assertEqual(len(result.columns), len(original_columns) + 1)

        # Verify values are positive
        self.assertTrue((result['area_test'] > 0).all())

    def test_get_km2_surface_area_without_inplace(self):
        """
        Test get_km2_surface_area method without inplace option.

        Verifies that with inplace=False, the method returns a new GeoDataFrame
        with an additional column containing area values in square kilometers.
        """
        result = self.gdf.get_km2_surface_area()
        self.assertIsInstance(result, StandardGeodataframe)
        self.assertIn('SURF_A_KM2', result.columns)

        # Verify values are positive
        self.assertTrue((result['SURF_A_KM2'] > 0).all())

    def test_get_km2_surface_area_with_inplace(self):
        """
        Test get_km2_surface_area method with inplace=True.

        Verifies that with inplace=True, the method adds a new column to
        the GeoDataFrame containing area values in square kilometers.
        """
        original_columns = self.gdf.columns.tolist()
        result = self.gdf.get_km2_surface_area(
            inplace=True, column_name='area_km2_test')
        self.assertIsInstance(result, StandardGeodataframe)
        self.assertIn('area_km2_test', result.columns)
        self.assertEqual(len(result.columns), len(original_columns) + 1)

        # Verify values are positive
        self.assertTrue((result['area_km2_test'] > 0).all())

    def test_calculate_surface_area_helper(self):
        """
        Test _calculate_surface_area helper method.

        Verifies that:
        1. The method returns a GeoDataFrame with a new area column
        2. The areas are calculated correctly according to the conversion factor
        3. The original CRS is preserved and returned
        """
        # First test with conversion factor of 1000 (m2 to km2)
        result, original_crs = self.gdf._calculate_surface_area(
            conversion_factor=1000, column_name='test_area')

        # Verify return type
        self.assertIsInstance(result, gpd.GeoDataFrame)

        # Verify that the column exists
        self.assertIn('test_area', result.columns)

        # Verify that original CRS is preserved
        self.assertEqual(original_crs, self.gdf.crs)

        # Verify that the area values are positive
        self.assertTrue((result['test_area'] > 0).all())

        # Second test with conversion factor of 1000000 (m2 to km2)
        result2, _ = self.gdf._calculate_surface_area(
            conversion_factor=1000000, column_name='test_area_km2')

        # Verify areas are proportional to the conversion factor
        for idx in range(len(result)):
            # Use assertAlmostEqual with a delta to allow rounding differences
            expected_ratio = result['test_area'].iloc[idx] / 1000
            actual_value = result2['test_area_km2'].iloc[idx]
            self.assertAlmostEqual(
                expected_ratio,
                actual_value,
                places=5,  # Allow rounding differences
                msg=f"Index failure {idx}: expected ~{expected_ratio}, obtained {actual_value}"
            )

    def test_substract_overlapping_geometries_with_tuple_args(self):
        """
        Test substract_overlapping_geometries with tuple arguments.

        Verifies that the method correctly processes tuple arguments
        in the format (target_value, values_to_subtract).
        """
        # Find a valid column for testing
        column_name = 'Name'  # Assuming this exists in the test data

        # Find valid values
        unique_values = self.gdf[column_name].unique()
        if len(unique_values) >= 2:
            target_value = unique_values[0]
            values_to_substract = [unique_values[1]]

            # Execute the function
            result = self.gdf.substract_overlapping_geometries(
                column_name,
                (target_value, values_to_substract)
            )
            self.assertIsInstance(result, StandardGeodataframe)

    def test_substract_overlapping_geometries_with_dict_args(self):
        """
        Test substract_overlapping_geometries with dictionary arguments.

        Verifies that the method correctly processes dictionary arguments
        in the format {target_value: [values_to_subtract]}.
        """
        # Find a valid column for testing
        column_name = 'Name'  # Assuming this exists in the test data

        # Find valid values
        unique_values = self.gdf[column_name].unique()
        if len(unique_values) >= 2:
            args = {unique_values[0]: [unique_values[1]]}

            # Execute the function
            result = self.gdf.substract_overlapping_geometries(
                column_name, args)
            self.assertIsInstance(result, StandardGeodataframe)

    def test_substract_overlapping_geometries_with_invalid_args(self):
        """
        Test substract_overlapping_geometries with invalid args format.

        Verifies that using an invalid arguments format raises ValueError.
        """
        column_name = 'Name'  # Assuming this exists in the test data

        # Invalid args format
        with self.assertRaises(ValueError):
            self.gdf.substract_overlapping_geometries(
                column_name, "invalid_args")

    def test_substract_overlapping_geometries_with_inplace(self):
        """
        Test substract_overlapping_geometries with inplace=True.

        Verifies that:
        1. With inplace=True, the original object is modified
        2. The same object instance is returned
        3. The result is a valid StandardGeodataframe
        """
        # Find a valid column for testing
        column_name = 'Name'  # Assuming this exists in the test data

        # Find valid values
        unique_values = self.gdf[column_name].unique()
        if len(unique_values) >= 2:
            args = {unique_values[0]: [unique_values[1]]}

            # Save reference to the original object
            original_id = id(self.gdf)

            # Execute the function with inplace=True
            result = self.gdf.substract_overlapping_geometries(
                column_name,
                args,
                inplace=True
            )

            # Verify that the the same instace is returned
            self.assertEqual(id(result), original_id)

            # Verify that result is the same object as self.gdf
            self.assertIs(result, self.gdf)

            # Verify that result is a StandardGeodataframe
            self.assertIsInstance(result, StandardGeodataframe)


class TestIntegration(unittest.TestCase):
    """
    Integration tests for the gdf_standardization module.

    This test class verifies the end-to-end functionality of the module,
    combining multiple operations from different classes to ensure they work
    together correctly in real-world scenarios.

    Attributes:
        test_file_path (str): Path to the test shapefile.
        gdf (StandardGeodataframe): A test StandardGeodataframe loaded from the shapefile.
        output_dir (str): Directory where output files will be saved for validation.
    """

    def setUp(self):
        """
        Set up test fixtures before each test method.

        Initializes:
        1. Path to a test shapefile
        2. A StandardGeodataframe loaded from the test file
        3. Output directory path for saving test results
        """
        current_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_file_path = os.path.join(
            current_directory,
            'tests_files',
            'inputs',
            'gdf_standardization_test_file.shp'
        )
        self.gdf = StandardGeodataframe(self.test_file_path)
        self.output_dir = os.path.join(
            current_directory,
            'tests_files',
            'outputs'
        )

    def test_end_to_end_workflow(self):
        """
        Test full workflow from file loading to geometry processing.

        Verifies that:
        1. A file can be loaded into a StandardGeodataframe
        2. CRS transformation is applied correctly
        3. Surface area calculation works on the transformed data
        4. Overlapping geometries can be subtracted from the result
        5. The final output maintains data integrity
        """
        # Step 1: Load data
        result = self.gdf
        self.assertIsInstance(result, StandardGeodataframe)

        # Step 2: Transform CRS
        result = result.to_crs("EPSG:3857", inplace=False)
        self.assertEqual(result.crs, "EPSG:3857")

        # Step 3: Calculate surface area
        result = StandardGeodataframe(result, crs="EPSG:3857")
        result = result.get_km2_surface_area(inplace=False)
        self.assertIn('SURF_A_KM2', result.columns)
        self.assertTrue(all(isinstance(x, (int, float))
                        for x in result['SURF_A_KM2']))

        # Step 4: Subtract overlapping geometries
        replacement_dict = {'3320844867': ['9041410429']}
        result = result.substract_overlapping_geometries(
            'Name', replacement_dict, inplace=False)

        # Step 5: Verify final output
        self.assertIsInstance(result, StandardGeodataframe)
        self.assertEqual(len(result), len(self.gdf))

    def test_geometry_validation_and_processing(self):
        """
        Test geometry validation and processing workflow.

        Verifies that:
        1. Invalid geometries are corrected during processing
        2. Z-coordinates are removed from 3D geometries
        3. Holes (geni) are removed from polygons
        4. Interior points can be calculated for the processed geometries
        """
        # Create a test GeoDataFrame with problematic geometries
        exterior = [(0, 0, 1), (0, 10, 1), (10, 10, 1), (10, 0, 1), (0, 0, 1)]
        hole = [(3, 3, 1), (3, 7, 1), (7, 7, 1), (7, 3, 1), (3, 3, 1)]
        with_hole = Polygon(exterior, [hole])

        # Self-intersecting polygon (invalid)
        invalid = Polygon([(0, 0), (10, 10), (0, 10), (10, 0), (0, 0)])

        test_gdf = gpd.GeoDataFrame(
            geometry=[with_hole, invalid], crs="EPSG:4326")

        # Process through StandardGeodataframe
        standard_gdf = StandardGeodataframe(test_gdf)
        result = standard_gdf

        # Verify Z-coordinates are removed
        for coord in result.geometry.iloc[0].exterior.coords:
            self.assertEqual(len(coord), 2)

        # Verify holes are removed
        self.assertEqual(len(result.geometry.iloc[0].interiors), 0)

        # Add interior points
        result = result.get_interior_points(inplace=True)
        self.assertIn('interior_point', result.columns)
        self.assertTrue(all(isinstance(pt, Point)
                        for pt in result['interior_point']))

    def test_import_export_roundtrip(self):
        """
        Test round-trip import/export functionality.

        Verifies that:
        1. A StandardGeodataframe can be saved to a file
        2. The file can be read back into a new StandardGeodataframe
        3. The data integrity is maintained through this process
        """
        # Skip if the test file doesn't exist
        if not os.path.exists(self.test_file_path):
            self.skipTest(f"Test file {self.test_file_path} not found.")

        # Step 1: Process original data
        original = self.gdf
        original_data = original.get_km2_surface_area(inplace=False)

        # Step 2: Save to file
        output_file = os.path.join(
            self.output_dir, "integration_test_roundtrip.shp")
        original_data.to_file(output_file)

        # Step 3: Read back and compare
        reimported = StandardGeodataframe(output_file)
        reimported_data = reimported

        # Compare structure and content
        self.assertEqual(set(original_data.columns),
                         set(reimported_data.columns))
        self.assertEqual(len(original_data), len(reimported_data))
        self.assertEqual(original_data.crs, reimported_data.crs)

        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)

    def test_combined_operations(self):
        """
        Test combining operations from different components.

        Verifies that:
        1. Input transformation works correctly with GeniRemover
        2. Adapter functionality integrates with StandardGeodataframe
        3. Multiple sequential operations produce consistent results
        """
        # Skip if the test file doesn't exist
        if not os.path.exists(self.test_file_path):
            self.skipTest(f"Test file {self.test_file_path} not found.")

        # Step 1: Create transformer and adapter
        transformer = _InputTransformer(self.test_file_path)
        adapter = GeoDataFrameAdapter(transformer.standard_gdf)

        # Step 2: Perform operations using the adapter
        result_gdf = adapter.to_crs("EPSG:3857")

        # Step 3: Create a new StandardGeodataframe from the result
        result = StandardGeodataframe(result_gdf, crs="EPSG:3857")

        # Step 4: Apply further processing
        result = result.get_km2_surface_area(inplace=False)
        result = result.get_interior_points(inplace=True)

        # Verify final outcome
        self.assertEqual(result.crs, "EPSG:3857")
        self.assertIn('SURF_A_KM2', result.columns)
        self.assertIn('interior_point', result.columns)


if __name__ == '__main__':
    unittest.main()
