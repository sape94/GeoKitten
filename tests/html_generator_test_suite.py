"""
This module contains unit tests and integration tests for the `html_generator` module.

It tests the functionality of the `_MapAttributeValidator`, `_CategoricalColorPaletteGenerator`, 
`_ContinuousColorPaletteGenerator`, `_MapStyleGenerator`, `_InteractiveMapCreator`, 
`InteractiveCategoricalHtmlMap`, and `InteractiveContinuousHtmlMap` classes.

The tests ensure proper handling of GeoDataFrames, color palette generation, map styling, 
and interactive map creation for both categorical and continuous data.

Author: Sergio A. Pelayo Escalera
Version: 0.1.0
Created: 2025-04-09
Last Updated: 2025-04-09
"""

from geokitten.html_generator import (
    _MapAttributeValidator,
    _CategoricalColorPaletteGenerator,
    _ContinuousColorPaletteGenerator,
    _MapStyleGenerator,
    _InteractiveMapCreator,
    InteractiveCategoricalHtmlMap,
    InteractiveContinuousHtmlMap
)
__version__ = "0.1.0"
__author__ = "Sergio A. Pelayo Escalera (sergio.pelayo@nielseniq.com)"
__collabs__ = [""]
__created__ = "2025-04-09"
__last_updated__ = "2025-04-09"

import os
import pytest
from unittest import mock
import geopandas as gpd
import folium
from shapely.geometry import Polygon
import branca.colormap as cm
import matplotlib.colors as mcolors
import sys
from pathlib import Path

# Add the src directory to the Python path
sys.path.append(str(Path(__file__).parent.parent))


# Test file paths
TEST_SHP_FILE = './tests/tests_files/inputs/html_generator_test_file.shp'
CATEGORICAL_OUTPUT = './tests/tests_files/outputs/html_generator_InteractiveCategoricalMap_test_validation_file.html'
CONTINUOUS_OUTPUT = './tests/tests_files/outputs/html_generator_InteractiveContinuousMap_test_validation_file.html'


def create_test_gdf():
    """
    Helper function to create a simple test GeoDataFrame.

    This function generates a GeoDataFrame with a single polygon and 
    associated attributes for testing purposes.

    Returns:
    --------
    GeoDataFrame
        A GeoDataFrame with one polygon and associated attributes.
    """
    polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    data = {
        'geometry': [polygon],
        'region': ['Region A'],
        'population': [1000],
        'id': [1],
        'name': ['Test Area']
    }
    return gpd.GeoDataFrame(data, crs="EPSG:4326")


def load_test_shapefile():
    """
    Helper function to load the test shapefile or create a mock GeoDataFrame.

    If the specified test shapefile exists, it is loaded as a GeoDataFrame.
    Otherwise, a synthetic GeoDataFrame is created with mock data.

    Returns:
    --------
    GeoDataFrame
        A GeoDataFrame loaded from the test shapefile or created with mock data.
    """
    if os.path.exists(TEST_SHP_FILE):
        return gpd.read_file(TEST_SHP_FILE)
    else:
        # If test file doesn't exist, create a mock GeoDataFrame
        print(
            f"Warning: Test file {TEST_SHP_FILE} not found. Using synthetic data.")
        gdf = create_test_gdf()
        # Add columns that would be in the real test file
        gdf['ID'] = gdf['id']
        gdf['DEPTO_ID'] = 1
        gdf['MUN_ID'] = 1
        gdf['MUN_NAME'] = 'Test Municipality'
        gdf['POPULATION'] = gdf['population']
        gdf['REGION'] = gdf['region']
        return gdf


class Test_MapAttributeValidator:
    """
    Unit tests for the `_MapAttributeValidator` class.

    These tests validate the functionality of GeoDataFrame validation, 
    color column validation, tooltip column validation, and map center 
    and bounds calculations.
    """

    def setup_method(self):
        """
        Set up test data before each test method.

        This method initializes a simple GeoDataFrame for testing.
        """
        self.gdf = create_test_gdf()

    def test_init_valid(self):
        """
        Test initialization with valid inputs.

        Ensures that the `_MapAttributeValidator` class initializes correctly 
        with a valid GeoDataFrame, color column, and tooltip columns.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region', 'population']
        )
        assert validator.gdf is self.gdf
        assert validator.color_column == 'region'
        assert validator.tooltip_columns == ['region', 'population']

    def test_init_invalid_gdf(self):
        """
        Test initialization with an invalid GeoDataFrame.

        Ensures that a `TypeError` is raised when the input is not a GeoDataFrame.
        """
        with pytest.raises(TypeError, match="Input must be a GeoDataFrame"):
            _MapAttributeValidator(
                {'data': 'not a geodataframe'}, 'region', ['region'])

    def test_init_invalid_color_column(self):
        """
        Test initialization with a non-existent color column.

        Ensures that a `ValueError` is raised when the specified color column 
        does not exist in the GeoDataFrame.
        """
        with pytest.raises(ValueError, match="Column 'nonexistent' not found in GeoDataFrame"):
            _MapAttributeValidator(self.gdf, 'nonexistent', ['region'])

    def test_init_invalid_tooltip_columns(self):
        """
        Test initialization with non-existent tooltip columns.

        Ensures that a `ValueError` is raised when none of the specified tooltip 
        columns exist in the GeoDataFrame.
        """
        with pytest.raises(ValueError, match="None of the specified tooltip columns exist in the GeoDataFrame"):
            _MapAttributeValidator(self.gdf, 'region', [
                                   'nonexistent1', 'nonexistent2'])

    def test_validate_geodataframe(self):
        """
        Test the `_validate_geodataframe` method.

        Ensures that the method correctly validates a GeoDataFrame and raises 
        a `TypeError` for invalid inputs.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region']
        )
        result = validator._validate_geodataframe(self.gdf)
        assert result is self.gdf

        with pytest.raises(TypeError, match="Input must be a GeoDataFrame"):
            validator._validate_geodataframe({'data': 'not a geodataframe'})

    def test_validate_color_column(self):
        """
        Test the `_validate_color_column` method.

        Ensures that the method correctly validates an existing color column 
        and raises a `ValueError` for non-existent columns.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region']
        )
        result = validator._validate_color_column('population')
        assert result == 'population'
        with pytest.raises(ValueError, match="Column 'nonexistent' not found in GeoDataFrame"):
            validator._validate_color_column('nonexistent')

    def test_validate_tooltip_columns(self):
        """
        Test the `_validate_tooltip_columns` method.

        Ensures that the method correctly validates and filters tooltip columns, 
        raising a `ValueError` if none of the specified columns exist.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region']
        )
        result = validator._validate_tooltip_columns(['population', 'id'])
        assert result == ['population', 'id']
        with pytest.raises(ValueError, match="None of the specified tooltip columns exist in the GeoDataFrame"):
            validator._validate_tooltip_columns(
                ['nonexistent1', 'nonexistent2'])

    def test_validate_tooltip_columns_warning(self, capsys):
        """
        Test warning when more than 8 tooltip columns are provided.

        Ensures that a warning is issued and only the first 8 columns are used.
        """
        # Create GDF with many columns
        data = {'geometry': [
            Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])]}
        cols = [f'col{i}' for i in range(10)]
        for i, col in enumerate(cols):
            data[col] = [i]
        gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
        validator = _MapAttributeValidator(gdf, 'col0', ['col0'])
        result = validator._validate_tooltip_columns(cols)
        # Check that only first 8 columns are returned
        assert len(result) == 8
        assert result == cols[:8]
        # Check warning output
        captured = capsys.readouterr()
        assert "Warning: Only the first 8 tooltip columns will be used" in captured.out

    def test_get_map_center_with_provided_location(self):
        """
        Test `_get_map_center` with a provided location.

        Ensures that the method returns the provided location without modification.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region']
        )
        location = [10.0, 20.0]
        result = validator._get_map_center(location)
        assert result == location

    def test_get_map_center_from_geodataframe(self):
        """
        Test `_get_map_center` with a calculated location from the GeoDataFrame.

        Ensures that the method calculates the centroid of the GeoDataFrame 
        when no location is provided.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region']
        )
        result = validator._get_map_center()
        # The centroid of the test polygon should be (0.5, 0.5)
        assert isinstance(result, list)
        assert len(result) == 2
        assert result[0] == 0.5  # latitude
        assert result[1] == 0.5  # longitude

    def test_prepare_geodata_already_geographic(self):
        """
        Test `_prepare_geodata` when the GeoDataFrame is already in geographic CRS.

        Ensures that the method returns the same GeoDataFrame without modification.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region']
        )
        result = validator._prepare_geodata()
        assert result is self.gdf  # Should return the same object

    def test_prepare_geodata_non_geographic(self):
        """
        Test `_prepare_geodata` when the GeoDataFrame is in a projected CRS.

        Ensures that the method converts the GeoDataFrame to a geographic CRS.
        """
        # Create a GeoDataFrame with a projected CRS
        gdf_projected = self.gdf.copy()
        gdf_projected.crs = "EPSG:3857"  # Web Mercator projection
        validator = _MapAttributeValidator(
            gdf_projected, 'region', ['region']
        )
        result = validator._prepare_geodata()
        # Should return a new GeoDataFrame with geographic CRS
        assert result is not gdf_projected
        assert result.crs.is_geographic

    def test_get_bounds(self):
        """
        Test the `_get_bounds` method.

        Ensures that the method correctly calculates the bounds of the GeoDataFrame 
        in geographic coordinates.
        """
        validator = _MapAttributeValidator(
            self.gdf, 'region', ['region']
        )
        bounds = validator._get_bounds()
        # Check structure and content of bounds
        assert isinstance(bounds, list)
        assert len(bounds) == 2
        assert len(bounds[0]) == 2
        assert len(bounds[1]) == 2
        # For our test polygon, bounds should be [[0, 0], [1, 1]]
        assert bounds[0][0] == 0  # min latitude
        assert bounds[0][1] == 0  # min longitude
        assert bounds[1][0] == 1  # max latitude
        assert bounds[1][1] == 1  # max longitude


class Test_CategoricalColorPaletteGenerator:
    """
    Unit tests for the `_CategoricalColorPaletteGenerator` class.

    These tests validate the functionality of generating color palettes 
    for categorical data, including custom colors, colormap usage, and 
    random color generation.
    """

    def setup_method(self):
        """
        Set up test data before each test method.

        This method initializes a `_CategoricalColorPaletteGenerator` instance 
        for testing.
        """
        self.generator = _CategoricalColorPaletteGenerator()

    def test_init_default(self):
        """
        Test initialization with default parameters.

        Ensures that the `_CategoricalColorPaletteGenerator` class initializes 
        with the default color scheme and no custom colors.
        """
        generator = _CategoricalColorPaletteGenerator()
        assert generator.color_scheme == 'tab20'
        assert generator.custom_colors is None

    def test_init_custom_scheme(self):
        """
        Test initialization with a custom color scheme.

        Ensures that the class initializes correctly with a specified color scheme.
        """
        generator = _CategoricalColorPaletteGenerator('Set1')
        assert generator.color_scheme == 'Set1'
        assert generator.custom_colors is None

    def test_validate_hex_colors_valid(self):
        """
        Test `_validate_hex_colors` with valid hex colors.

        Ensures that the method correctly validates a list of valid hex color codes.
        """
        colors = ['#FF0000', '#00FF00', '#0000FF']
        result = self.generator._validate_hex_colors(colors)
        assert result == colors

    def test_validate_hex_colors_invalid(self):
        """
        Test `_validate_hex_colors` with invalid hex colors.

        Ensures that the method raises a `ValueError` for invalid hex color codes.
        """
        # Missing # prefix
        with pytest.raises(ValueError, match="Invalid hex color code: FF0000"):
            self.generator._validate_hex_colors(['FF0000'])
        # Invalid characters
        with pytest.raises(ValueError, match="Invalid hex color code: #XYZ123"):
            self.generator._validate_hex_colors(['#XYZ123'])

    def test_set_custom_colors(self):
        """
        Test `_set_custom_colors` method.

        Ensures that the method correctly sets custom colors and updates the 
        color scheme to 'custom'.
        """
        colors = ['#FF0000', '#00FF00', '#0000FF']
        self.generator._set_custom_colors(colors)
        assert self.generator.custom_colors == colors
        assert self.generator.color_scheme == 'custom'

    def test_set_custom_colors_none(self):
        """
        Test `_set_custom_colors` with `None`.

        Ensures that the method resets the custom colors and retains the default 
        color scheme.
        """
        self.generator._set_custom_colors(None)
        assert self.generator.custom_colors is None
        assert self.generator.color_scheme == 'tab20'  # Should keep default

    def test_generate_distinct_colors(self):
        """
        Test `_generate_distinct_colors` method.

        Ensures that the method generates the specified number of distinct colors 
        and all colors are valid hex codes.
        """
        # Generate 5 distinct colors
        colors = self.generator._generate_distinct_colors(5)
        assert len(colors) == 5
        # All colors should be unique
        assert len(set(colors)) == 5
        # All colors should be valid hex codes
        for color in colors:
            assert color.startswith('#')
            assert len(color) == 7

    def test_use_custom_colors_enough(self):
        """
        Test `_use_custom_colors` when enough custom colors are provided.

        Ensures that the method returns the specified number of custom colors.
        """
        self.generator.custom_colors = ['#FF0000', '#00FF00', '#0000FF']
        result = self.generator._use_custom_colors(2)
        assert result == ['#FF0000', '#00FF00']

    def test_use_custom_colors_not_enough(self):
        """
        Test `_use_custom_colors` when not enough custom colors are provided.

        Ensures that the method cycles through the custom colors to meet the 
        required number.
        """
        self.generator.custom_colors = ['#FF0000', '#00FF00']
        result = self.generator._use_custom_colors(5)
        assert result == ['#FF0000', '#00FF00',
                          '#FF0000', '#00FF00', '#FF0000']

    def test_try_specific_colormap_valid(self):
        """
        Test `_try_specific_colormap` with a valid colormap.

        Ensures that the method returns a list of colors from the specified colormap.
        """
        self.generator.color_scheme = 'viridis'
        result = self.generator._try_specific_colormap(3)
        assert len(result) == 3
        assert all(color.startswith('#') for color in result)

    def test_try_specific_colormap_invalid(self):
        """
        Test `_try_specific_colormap` with an invalid colormap.

        Ensures that the method returns `None` for an invalid colormap.
        """
        self.generator.color_scheme = 'nonexistent_colormap'
        result = self.generator._try_specific_colormap(3)
        assert result is None

    def test_collect_colors_from_multiple_sources(self):
        """
        Test `_collect_colors_from_multiple_sources` method.

        Ensures that the method collects a sufficient number of colors from 
        various sources and all colors are valid hex codes.
        """
        result = self.generator._collect_colors_from_multiple_sources(50)
        # Should return a significant number of colors
        assert len(result) > 20
        # All should be valid hex codes
        assert all(color.startswith('#') for color in result)

    def test_add_colors_from_colormaps(self):
        """
        Test `_add_colors_from_colormaps` method.

        Ensures that the method adds colors from various colormaps to the list.
        """
        colors = []
        self.generator._add_colors_from_colormaps(colors, 10)
        # Should have added colors from various colormaps
        assert len(colors) > 0
        # All should be valid hex codes
        assert all(color.startswith('#') for color in colors)

    def test_add_random_colors(self):
        """
        Test `_add_random_colors` method.

        Ensures that the method generates random colors and appends them to the 
        existing list of colors.
        """
        existing_colors = ['#FF0000', '#00FF00']
        result = self.generator._add_random_colors(existing_colors, 5)
        # Should have 5 colors in total
        assert len(result) == 5
        # Original colors should be preserved
        assert result[:2] == existing_colors
        # All should be valid hex codes
        assert all(color.startswith('#') for color in result)

    def test_generate_distinct_random_color(self):
        """
        Test `_generate_distinct_random_color` method.

        Ensures that the method generates a valid hex color that is distinct 
        from existing colors.
        """
        existing_hsv = []  # Empty list, so any color should be distinct
        result = self.generator._generate_distinct_random_color(existing_hsv)
        # Should return a valid hex color
        assert result.startswith('#')
        assert len(result) == 7

    def test_is_color_distinct_true(self):
        """
        Test `_is_color_distinct` when the color is distinct.

        Ensures that the method returns `True` for a color that is sufficiently 
        distinct from existing colors.
        """
        # Create HSV values representing red
        existing_hsv = [[0, 1, 1]]  # Red in HSV
        # Test with blue (h=0.66, s=1, v=1), which should be distinct from red
        result = self.generator._is_color_distinct(0.66, 1, 1, existing_hsv)
        assert result is True

    def test_is_color_distinct_false(self):
        """
        Test `_is_color_distinct` when the color is not distinct.

        Ensures that the method returns `False` for a color that is not sufficiently 
        distinct from existing colors.
        """
        # Create HSV values representing red
        existing_hsv = [[0, 1, 1]]  # Red in HSV
        # Test with slightly different red, which should not be distinct enough
        result = self.generator._is_color_distinct(
            0.01, 0.99, 0.99, existing_hsv)
        assert result is False

    def test_create_categorical_color_dict_numeric(self):
        """
        Test `_create_categorical_color_dict` with numeric values.

        Ensures that the method creates a dictionary mapping numeric values to 
        distinct colors.
        """
        values = [3, 1, 2]
        result = self.generator._create_categorical_color_dict(values)
        # Should return a dictionary with sorted keys
        assert set(result.keys()) == set(values)
        assert list(result.keys()) == [1, 2, 3]  # Should be sorted
        # Values should be valid hex colors
        assert all(color.startswith('#') for color in result.values())

    def test_create_categorical_color_dict_strings(self):
        """
        Test `_create_categorical_color_dict` with string values.

        Ensures that the method creates a dictionary mapping string values to 
        distinct colors.
        """
        values = ['c', 'a', 'b']
        result = self.generator._create_categorical_color_dict(values)
        # For strings, order in the result should match order of unique values in input
        assert set(result.keys()) == set(values)
        # Values should be valid hex colors
        assert all(color.startswith('#') for color in result.values())

    def test_create_categorical_color_dict_single_custom_color(self):
        """
        Test `_create_categorical_color_dict` with a single custom color.

        Ensures that all values are mapped to the same custom color.
        """
        self.generator.custom_colors = ['#FF0000']
        self.generator.color_scheme = 'custom'
        values = ['a', 'b', 'c']
        result = self.generator._create_categorical_color_dict(values)
        # All values should map to the same color
        assert set(result.values()) == {'#FF0000'}
        assert result == {'a': '#FF0000', 'b': '#FF0000', 'c': '#FF0000'}


class Test_ContinuousColorPaletteGenerator:
    """
    Unit tests for the `_ContinuousColorPaletteGenerator` class.

    These tests validate the functionality of generating continuous 
    colormaps for numeric data, including colormap validation and 
    color range handling.
    """

    def setup_method(self):
        """
        Set up test data before each test method.

        This method initializes a `_ContinuousColorPaletteGenerator` instance 
        with default parameters for testing.
        """
        self.generator = _ContinuousColorPaletteGenerator()

    def test_init_default(self):
        """
        Test initialization with default parameters.

        Ensures that the `_ContinuousColorPaletteGenerator` class initializes 
        with the default color scheme ('viridis').
        """
        generator = _ContinuousColorPaletteGenerator()
        assert generator.color_scheme == 'viridis'

    def test_init_custom_scheme(self):
        """
        Test initialization with a custom color scheme.

        Ensures that the class initializes correctly with a specified valid 
        color scheme (e.g., 'plasma').
        """
        generator = _ContinuousColorPaletteGenerator('plasma')
        assert generator.color_scheme == 'plasma'

    def test_init_invalid_scheme(self):
        """
        Test initialization with an invalid color scheme.

        Ensures that a `ValueError` is raised when an invalid color scheme 
        is provided during initialization.
        """
        with pytest.raises(ValueError, match="Color scheme 'nonexistent' not found in matplotlib"):
            _ContinuousColorPaletteGenerator('nonexistent')

    def test_validate_color_scheme_valid(self):
        """
        Test `_validate_color_scheme` with a valid color scheme.

        Ensures that the method validates a valid color scheme without 
        raising any exceptions.
        """
        # This should not raise an exception
        self.generator._validate_color_scheme()

    def test_validate_color_scheme_invalid(self):
        """Test _validate_color_scheme with invalid color scheme."""
        self.generator.color_scheme = 'nonexistent'
        with pytest.raises(ValueError, match="Color scheme 'nonexistent' not found in matplotlib"):
            self.generator._validate_color_scheme()

    def test_create_continuous_colormap(self):
        """
        Test `_validate_color_scheme` with an invalid color scheme.

        Ensures that a `ValueError` is raised when the color scheme is not 
        found in Matplotlib's colormaps.
        """
        min_value = 0
        max_value = 100
        result = self.generator._create_continuous_colormap(
            min_value, max_value)
        # Should return a LinearColormap
        assert isinstance(result, cm.LinearColormap)
        # Check colormap properties
        assert result.vmin == min_value
        assert result.vmax == max_value
        # Should have 256 colors
        assert len(result.colors) == 256
        # Colors should be valid hex codes
        assert all(mcolors.to_hex(color).startswith('#')
                   for color in result.colors)


class Test_MapStyleGenerator:
    """
    Unit tests for the `_MapStyleGenerator` class.

    These tests validate the functionality of creating style functions 
    and legends for both categorical and continuous data.
    """

    def setup_method(self):
        """
        Set up test data before each test method.

        This method initializes a `_MapStyleGenerator` instance with default 
        parameters for testing.
        """
        self.generator = _MapStyleGenerator()

    def test_init_default(self):
        """
        Test initialization with default parameters.

        Ensures that the `_MapStyleGenerator` class initializes with the default 
        values for `fill_opacity`, `line_opacity`, and `line_weight`.
        """
        generator = _MapStyleGenerator()
        assert generator.fill_opacity == 0.5
        assert generator.line_opacity == 1.0
        assert generator.line_weight == 1

    def test_init_custom(self):
        """
        Test initialization with custom parameters.

        Ensures that the `_MapStyleGenerator` class initializes correctly with 
        custom values for `fill_opacity`, `line_opacity`, and `line_weight`.
        """
        generator = _MapStyleGenerator(
            fill_opacity=0.7, line_opacity=0.8, line_weight=2)
        assert generator.fill_opacity == 0.7
        assert generator.line_opacity == 0.8
        assert generator.line_weight == 2

    def test_create_categorical_style(self):
        """
        Test `_create_categorical_style` method.

        Ensures that the method generates a valid style function and legend HTML 
        for categorical data, using a provided color column and color dictionary.
        """
        color_column = 'region'
        color_dict = {'Region A': '#FF0000', 'Region B': '#00FF00'}
        style_function, legend_html = self.generator._create_categorical_style(
            color_column, color_dict)
        # Test style function
        feature = {'properties': {'region': 'Region A'}}
        style = style_function(feature)
        assert style['fillColor'] == '#FF0000'
        assert style['color'] == 'black'
        assert style['weight'] == 1
        assert style['fillOpacity'] == 0.5
        assert style['opacity'] == 1.0
        # Test legend HTML
        assert isinstance(legend_html, str)
        assert 'region' in legend_html
        assert '#FF0000' in legend_html
        assert 'Region A' in legend_html
        assert '#00FF00' in legend_html
        assert 'Region B' in legend_html

    def test_create_categorical_style_with_missing_value(self):
        """
        Test `_create_categorical_style` with a feature that has a missing value.

        Ensures that the style function assigns a default white color (`#FFFFFF`) 
        to features with values not present in the color dictionary.
        """
        color_column = 'region'
        color_dict = {'Region A': '#FF0000', 'Region B': '#00FF00'}
        style_function, _ = self.generator._create_categorical_style(
            color_column, color_dict)
        # Test with a feature that has a region not in the color_dict
        feature = {'properties': {'region': 'Region C'}}
        style = style_function(feature)
        assert style['fillColor'] == '#FFFFFF'  # Default white color

    def test_create_continuous_style(self):
        """
        Test `_create_continuous_style` method.

        Ensures that the method generates a valid style function for continuous 
        data, using a provided color column and colormap.
        """
        color_column = 'population'
        # Create a mock colormap
        colormap = mock.MagicMock()
        colormap.return_value = '#FF0000'  # Always return red
        style_function = self.generator._create_continuous_style(
            color_column, colormap)
        # Test style function
        feature = {'properties': {'population': 1000}}
        style = style_function(feature)
        assert style['fillColor'] == '#FF0000'
        assert style['color'] == 'black'
        assert style['weight'] == 1
        assert style['fillOpacity'] == 0.5
        assert style['opacity'] == 1.0
        # Verify colormap was called with the correct value
        colormap.assert_called_once_with(1000)


class Test_InteractiveMapCreator:
    """
    Unit tests for the `_InteractiveMapCreator` class.

    These tests validate the functionality of creating base maps, 
    adding choropleth layers, legends, controls, and titles to the map.
    """

    def setup_method(self):
        """
        Set up test data before each test method.

        This method initializes a `_MapAttributeValidator` instance and an 
        `_InteractiveMapCreator` instance with a test GeoDataFrame for testing.
        """
        self.gdf = create_test_gdf()
        self.validator = _MapAttributeValidator(
            self.gdf, 'region', ['region', 'population'])
        self.map_creator = _InteractiveMapCreator(self.validator)

    def test_init(self):
        """
        Test initialization of the `_InteractiveMapCreator` class.

        Ensures that the class initializes correctly with a valid `_MapAttributeValidator` instance.
        """
        map_creator = _InteractiveMapCreator(self.validator)
        assert map_creator.validator is self.validator

    def test_create_map_default(self):
        """
        Test the `_create_map` method with default parameters.

        Ensures that the method creates a folium map with auto-calculated bounds 
        and default settings.
        """
        m = self.map_creator._create_map()
        assert isinstance(m, folium.Map)
        # Should auto-fit to bounds
        assert m.get_bounds() is not None

    def test_create_map_with_location_and_zoom(self):
        """
        Test the `_create_map` method with specified location and zoom level.

        Ensures that the method creates a folium map with the provided center 
        location and zoom level.
        """
        location = [10.0, 20.0]
        zoom_start = 8
        m = self.map_creator._create_map(
            zoom_start=zoom_start, location=location)
        assert isinstance(m, folium.Map)
        # Map should have the specified location and zoom
        assert m.location == location
        assert m.options['zoom'] == zoom_start

    def test_add_choropleth_layer(self):
        """
        Test the `_add_choropleth_layer` method.

        Ensures that the method adds a GeoJSON layer with the specified style 
        function and tooltip columns to the folium map.
        """
        m = folium.Map()
        def style_function(x): return {'fillColor': '#FF0000'}
        tooltip_columns = ['region', 'population']
        result = self.map_creator._add_choropleth_layer(
            m, style_function, tooltip_columns)
        assert result is m  # Should return the same map object
        # Check that GeoJson layer was added
        assert len(m._children) > 0
        assert any(isinstance(child, folium.GeoJson)
                   for child in m._children.values())

    def test_add_controls(self):
        """
        Test the `_add_controls` method.

        Ensures that the method adds a `LayerControl` to the folium map.
        """
        m = folium.Map()
        result = self.map_creator._add_controls(m)
        assert result is m  # Should return the same map object
        # Check that LayerControl was added
        assert any(isinstance(child, folium.LayerControl)
                   for child in m._children.values())

    def test_add_title(self):
        """
        Test the `_add_title` method.

        Ensures that the method adds a title to the folium map's HTML structure.
        """
        m = folium.Map()
        title = "Test Map"
        result = self.map_creator._add_title(m, title)
        assert result is m  # Should return the same map object
        # Title should be added to HTML
        html = m.get_root().render()
        assert title in html

    def test_add_title_none(self):
        """
        Test the `_add_title` method with a `None` title.

        Ensures that the method does not add a title to the folium map when 
        no title is provided.
        """
        m = folium.Map()
        result = self.map_creator._add_title(m, None)
        assert result is m  # Should return the same map object
        # No title element should be added
        html = m.get_root().render()
        assert "<h3" not in html

    def test_add_categorical_legend(self):
        """
        Test the `_add_categorical_legend` method.

        Ensures that the method adds a categorical legend to the folium map's 
        HTML structure.
        """
        m = folium.Map()
        legend_html = "<div>Test Legend</div>"
        result = self.map_creator._add_categorical_legend(m, legend_html)
        assert result is m  # Should return the same map object
        # Legend should be added to HTML
        html = m.get_root().render()
        assert legend_html in html

    def test_add_continuous_colormap(self):
        """
        Test the `_add_continuous_colormap` method.

        Ensures that the method adds a continuous colormap to the folium map 
        and sets the colormap's caption.
        """
        m = folium.Map()
        colormap = cm.LinearColormap(['red', 'blue'], vmin=0, vmax=100)
        legend_name = "Test Colormap"
        result = self.map_creator._add_continuous_colormap(
            m, colormap, legend_name)
        assert result is m  # Should return the same map object
        # Colormap should be added to the map
        assert any(isinstance(child, cm.LinearColormap)
                   for child in m._children.values())
        # Caption should be set
        for child in m._children.values():
            if isinstance(child, cm.LinearColormap):
                assert child.caption == legend_name


class TestInteractiveCategoricalHtmlMap:
    """
    Unit tests for the `InteractiveCategoricalHtmlMap` class.

    These tests validate the functionality of creating categorical 
    interactive maps, including color palette generation, map styling, 
    and legend creation.
    """

    def setup_method(self):
        """
        Set up test data before each test method.

        This method initializes a test GeoDataFrame and an `InteractiveCategoricalHtmlMap` 
        instance for testing.
        """
        self.gdf = create_test_gdf()
        self.map_creator = InteractiveCategoricalHtmlMap(
            self.gdf, 'region', ['region', 'population']
        )

    def test_init(self):
        """
        Test initialization of the `InteractiveCategoricalHtmlMap` class.

        Ensures that the class initializes correctly with a valid GeoDataFrame, 
        color column, and tooltip columns.
        """
        map_creator = InteractiveCategoricalHtmlMap(
            self.gdf, 'region', ['region', 'population']
        )
        assert isinstance(map_creator.validator, _MapAttributeValidator)
        assert isinstance(map_creator.map_creator, _InteractiveMapCreator)

    def test_initialize_generators(self):
        """
        Test the `_initialize_generators` method.

        Ensures that the method initializes the color palette generator and 
        style generator with the specified parameters.
        """
        color_generator, style_generator = self.map_creator._initialize_generators(
            color_scheme='Set1',
            set_custom_colors=['#FF0000', '#00FF00'],
            fill_opacity=0.7,
            line_opacity=0.8,
            line_weight=2
        )
        assert isinstance(color_generator, _CategoricalColorPaletteGenerator)
        assert color_generator.color_scheme == 'custom'  # Changed by set_custom_colors
        assert color_generator.custom_colors == ['#FF0000', '#00FF00']
        assert isinstance(style_generator, _MapStyleGenerator)
        assert style_generator.fill_opacity == 0.7
        assert style_generator.line_opacity == 0.8
        assert style_generator.line_weight == 2

    def test_create_map(self):
        """
        Test the `_create_map` method.

        Ensures that the method calls the `_create_map` method of the 
        `_InteractiveMapCreator` class with the correct parameters.
        """
        location = [10.0, 20.0]
        zoom_start = 8
        # Mock the _create_map method of _InteractiveMapCreator
        with mock.patch.object(self.map_creator.map_creator, '_create_map') as mock_create_map:
            self.map_creator._create_map(location, zoom_start)
            mock_create_map.assert_called_once_with(zoom_start, location)

    def test_apply_categorical_styling(self):
        """
        Test the `_apply_categorical_styling` method.

        Ensures that the method generates a valid style function and legend HTML 
        for categorical data using the specified color generator and style generator.
        """
        color_generator = _CategoricalColorPaletteGenerator()
        style_generator = _MapStyleGenerator()
        legend_name = "Test Legend"
        style_function, legend_html = self.map_creator._apply_categorical_styling(
            color_generator, style_generator, legend_name
        )
        assert callable(style_function)
        assert isinstance(legend_html, str)

    def test_add_layers_controls(self):
        """
        Test the `_add_layers_controls` method.

        Ensures that the method adds layers, controls, and a title to the map 
        using the specified style function, legend HTML, and title.
        """
        m = folium.Map()
        def style_function(x): return {'fillColor': '#FF0000'}
        legend_html = "<div>Test Legend</div>"
        title = "Test Map"
        # Mock the methods of _InteractiveMapCreator
        with mock.patch.object(self.map_creator.map_creator, '_add_choropleth_layer') as mock_add_layer, \
                mock.patch.object(self.map_creator.map_creator, '_add_categorical_legend') as mock_add_legend, \
                mock.patch.object(self.map_creator.map_creator, '_add_controls') as mock_add_controls, \
                mock.patch.object(self.map_creator.map_creator, '_add_title') as mock_add_title:

            mock_add_layer.return_value = m
            mock_add_legend.return_value = m
            mock_add_controls.return_value = m
            mock_add_title.return_value = m
            result = self.map_creator._add_layers_controls(
                m, style_function, legend_html, title)
            assert result is m
            mock_add_layer.assert_called_once()
            mock_add_legend.assert_called_once_with(m, legend_html)
            mock_add_controls.assert_called_once()
            mock_add_title.assert_called_once_with(m, title)

    def test_create(self, monkeypatch):
        """
        Test the `create` method.

        Ensures that the method creates a categorical interactive map with the 
        specified parameters and saves it to the specified output file.
        """
        # Mock the underlying methods
        mock_map = mock.MagicMock(spec=folium.Map)
        monkeypatch.setattr(self.map_creator, '_initialize_generators',
                            lambda *args, **kwargs: (mock.MagicMock(), mock.MagicMock()))
        monkeypatch.setattr(self.map_creator, '_create_map',
                            lambda *args, **kwargs: mock_map)
        monkeypatch.setattr(self.map_creator, '_apply_categorical_styling',
                            lambda *args, **kwargs: (lambda x: None, "legend"))
        monkeypatch.setattr(self.map_creator, '_add_layers_controls',
                            lambda *args, **kwargs: mock_map)
        # Mock the save method of the map
        mock_map.save = mock.MagicMock()
        result = self.map_creator.create(
            color_scheme='Set1',
            set_custom_colors=['#FF0000', '#00FF00'],
            fill_opacity=0.7,
            line_opacity=0.8,
            line_weight=2,
            location=[10.0, 20.0],
            zoom_start=8,
            title="Test Map",
            output_file="html_test_temp_dir/test_map.html"
        )
        assert result is mock_map
        mock_map.save.assert_called_once_with(
            "html_test_temp_dir/test_map.html")
        os.rmdir("html_test_temp_dir")

    def test_integration_with_real_data(self):
        """
        Integration test with real-world data.

        Ensures that the `InteractiveCategoricalHtmlMap` class can create a 
        categorical map using a real-world dataset and specified parameters.
        """
        # Skip if shapefile doesn't exist
        gdf = load_test_shapefile()
        map_creator = InteractiveCategoricalHtmlMap(
            gdf=gdf,
            color_column='REGION',
            tooltip_columns=['ID', 'DEPTO_ID', 'MUN_ID',
                             'MUN_NAME', 'POPULATION', 'REGION']
        )
        # Use MagicMock to avoid actually writing to disk
        with mock.patch.object(folium.Map, 'save'):
            map_obj = map_creator.create(
                fill_opacity=0.6,
                set_custom_colors=None,
                title='Medellín - Regions',
                output_file=CATEGORICAL_OUTPUT
            )
            # Verify the result is a folium Map
            assert isinstance(map_obj, folium.Map)
            # Verify it has GeoJson layers
            assert any(isinstance(child, folium.GeoJson)
                       for child in map_obj._children.values())


class TestInteractiveContinuousHtmlMap:
    """
    Unit tests for the `InteractiveContinuousHtmlMap` class.

    These tests validate the functionality of creating continuous 
    interactive maps, including colormap generation, map styling, 
    and legend creation.
    """

    def setup_method(self):
        """
        Set up test data before each test method.

        This method initializes a test GeoDataFrame and an `InteractiveContinuousHtmlMap` 
        instance for testing.
        """
        self.gdf = create_test_gdf()
        self.map_creator = InteractiveContinuousHtmlMap(
            self.gdf, 'population', ['region', 'population']
        )

    def test_init(self):
        """
        Test initialization of the `InteractiveContinuousHtmlMap` class.

        Ensures that the class initializes correctly with a valid GeoDataFrame, 
        numeric color column, and tooltip columns.
        """
        map_creator = InteractiveContinuousHtmlMap(
            self.gdf, 'population', ['region', 'population']
        )
        assert isinstance(map_creator.validator, _MapAttributeValidator)
        assert isinstance(map_creator.map_creator, _InteractiveMapCreator)

    def test_init_with_non_numeric_column(self):
        """
        Test initialization with a non-numeric color column.

        Ensures that a `ValueError` is raised when the specified color column 
        is not numeric.
        """
        with pytest.raises(ValueError, match="Color column 'region' must be numeric for continuous maps"):
            InteractiveContinuousHtmlMap(
                self.gdf, 'region', ['region', 'population']
            )

    def test_validate_color_column(self):
        """
        Test the `_validate_color_column` method.

        Ensures that the method validates a numeric color column without raising 
        an exception and raises a `ValueError` for non-numeric columns.
        """
        # Should not raise an exception for numeric column
        self.map_creator._validate_color_column(self.gdf, 'population')
        # Should raise an exception for non-numeric column
        with pytest.raises(ValueError, match="Color column 'region' must be numeric for continuous maps"):
            self.map_creator._validate_color_column(self.gdf, 'region')

    def test_initialize_generators(self):
        """
        Test the `_initialize_generators` method.

        Ensures that the method initializes the continuous color palette generator 
        and style generator with the specified parameters.
        """
        color_generator, style_generator = self.map_creator._initialize_generators(
            color_scheme='plasma',
            fill_opacity=0.7,
            line_opacity=0.8,
            line_weight=2
        )
        assert isinstance(color_generator, _ContinuousColorPaletteGenerator)
        assert color_generator.color_scheme == 'plasma'
        assert isinstance(style_generator, _MapStyleGenerator)
        assert style_generator.fill_opacity == 0.7
        assert style_generator.line_opacity == 0.8
        assert style_generator.line_weight == 2

    def test_create_map(self):
        """
        Test the `_create_map` method.

        Ensures that the method calls the `_create_map` method of the 
        `_InteractiveMapCreator` class with the correct parameters.
        """
        location = [10.0, 20.0]
        zoom_start = 8
        # Mock the _create_map method of _InteractiveMapCreator
        with mock.patch.object(self.map_creator.map_creator, '_create_map') as mock_create_map:
            self.map_creator._create_map(location, zoom_start)
            mock_create_map.assert_called_once_with(zoom_start, location)

    def test_apply_continuous_styling(self):
        """
        Test the `_apply_continuous_styling` method.

        Ensures that the method generates a valid style function and colormap 
        for continuous data using the specified color generator and style generator.
        """
        color_generator = _ContinuousColorPaletteGenerator()
        style_generator = _MapStyleGenerator()
        legend_name = "Test Legend"
        style_function, colormap = self.map_creator._apply_continuous_styling(
            color_generator, style_generator, legend_name
        )
        assert callable(style_function)
        assert isinstance(colormap, cm.LinearColormap)
        # Caption is set in _add_continuous_colormap, not here
        assert colormap.caption != legend_name

    def test_add_layers_controls(self):
        """
        Test the `_add_layers_controls` method.

        Ensures that the method adds layers, controls, and a title to the map 
        using the specified style function, colormap, legend name, and title.
        """
        m = folium.Map()
        def style_function(x): return {'fillColor': '#FF0000'}
        colormap = cm.LinearColormap(['red', 'blue'], vmin=0, vmax=100)
        legend_name = "Test Legend"
        title = "Test Map"
        # Mock the methods of _InteractiveMapCreator
        with mock.patch.object(self.map_creator.map_creator, '_add_choropleth_layer') as mock_add_layer, \
                mock.patch.object(self.map_creator.map_creator, '_add_continuous_colormap') as mock_add_colormap, \
                mock.patch.object(self.map_creator.map_creator, '_add_controls') as mock_add_controls, \
                mock.patch.object(self.map_creator.map_creator, '_add_title') as mock_add_title:

            mock_add_layer.return_value = m
            mock_add_colormap.return_value = m
            mock_add_controls.return_value = m
            mock_add_title.return_value = m
            result = self.map_creator._add_layers_controls(
                m, style_function, colormap, legend_name, title)
            assert result is m
            mock_add_layer.assert_called_once()
            mock_add_colormap.assert_called_once_with(m, colormap, legend_name)
            mock_add_controls.assert_called_once()
            mock_add_title.assert_called_once_with(m, title)

    def test_create(self, monkeypatch):
        """
        Test the `create` method.

        Ensures that the method creates a continuous interactive map with the 
        specified parameters and saves it to the specified output file.
        """
        # Mock the underlying methods
        mock_map = mock.MagicMock(spec=folium.Map)
        monkeypatch.setattr(self.map_creator, '_initialize_generators',
                            lambda *args, **kwargs: (mock.MagicMock(), mock.MagicMock()))
        monkeypatch.setattr(self.map_creator, '_create_map',
                            lambda *args, **kwargs: mock_map)
        monkeypatch.setattr(self.map_creator, '_apply_continuous_styling',
                            lambda *args, **kwargs: (lambda x: None, mock.MagicMock()))
        monkeypatch.setattr(self.map_creator, '_add_layers_controls',
                            lambda *args, **kwargs: mock_map)
        # Mock the save method of the map
        mock_map.save = mock.MagicMock()
        result = self.map_creator.create(
            color_scheme='plasma',
            fill_opacity=0.7,
            line_opacity=0.8,
            line_weight=2,
            location=[10.0, 20.0],
            zoom_start=8,
            title="Test Map",
            legend_name="Test Legend",
            output_file="html_test_temp_dir/test_map.html"
        )
        assert result is mock_map
        mock_map.save.assert_called_once_with(
            "html_test_temp_dir/test_map.html")
        os.rmdir("html_test_temp_dir")

    def test_integration_with_real_data(self):
        """
        Integration test with real-world data.

        Ensures that the `InteractiveContinuousHtmlMap` class can create a 
        continuous map using a real-world dataset and specified parameters.
        """
        # Skip if shapefile doesn't exist
        gdf = load_test_shapefile()
        map_creator = InteractiveContinuousHtmlMap(
            gdf=gdf,
            color_column='POPULATION',
            tooltip_columns=['ID', 'DEPTO_ID', 'MUN_ID',
                             'MUN_NAME', 'POPULATION', 'REGION']
        )
        # Use MagicMock to avoid actually writing to disk
        with mock.patch.object(folium.Map, 'save'):
            map_obj = map_creator.create(
                fill_opacity=0.6,
                title='Medellín - Population',
                output_file=CONTINUOUS_OUTPUT
            )
            # Verify the result is a folium Map
            assert isinstance(map_obj, folium.Map)
            # Verify it has GeoJson layers
            assert any(isinstance(child, folium.GeoJson)
                       for child in map_obj._children.values())
            # Verify it has a ColorMap
            assert any(isinstance(child, cm.LinearColormap)
                       for child in map_obj._children.values())


class TestIntegration:
    """
    Integration tests for the `html_generator` module.

    These tests validate the full workflow of the `InteractiveCategoricalHtmlMap` 
    and `InteractiveContinuousHtmlMap` classes with real-world data and 
    specific visualization parameters.
    """

    @pytest.mark.skipif(not os.path.exists(os.path.join("tests", "tests_files", "inputs", "html_generator_test_file.shp")),
                        reason="Test file html_generator_test_file.shp not found")
    def test_categorical_map_workflow(self):
        """
        Test the full workflow for categorical map creation.

        Ensures that the categorical map is correctly created with the specified
        parameters and can be rendered without errors.
        """
        # Define input and output file paths
        input_file = os.path.join(
            "tests", "tests_files", "inputs", "html_generator_test_file.shp")
        output_file = os.path.join("tests", "tests_files", "outputs",
                                   "html_generator_InteractiveCategoricalMap_test_validation_file.html")
        # Load the test file
        gdf = gpd.read_file(input_file)
        # Initialize map creator with appropriate columns
        map_creator = InteractiveCategoricalHtmlMap(
            gdf=gdf,
            color_column='REGION',
            tooltip_columns=['ID', 'DEPTO_ID', 'MUN_ID',
                             'MUN_NAME', 'POPULATION', 'REGION']
        )
        # Create the map without saving it
        with mock.patch('folium.Map.save'):
            map_obj = map_creator.create(
                fill_opacity=0.6,
                set_custom_colors=None,
                title='Medellín - Regions',
                output_file=output_file
            )
        # Verify the map was created successfully
        assert isinstance(map_obj, folium.Map)
        # Verify it has the expected components
        assert any(isinstance(child, folium.features.GeoJson)
                   for child in map_obj._children.values())
        # Verify tooltip columns are present
        geo_json_layers = [child for child in map_obj._children.values()
                           if isinstance(child, folium.features.GeoJson)]

    @pytest.mark.skipif(not os.path.exists(os.path.join("tests", "tests_files", "inputs", "html_generator_test_file.shp")),
                        reason="Test file html_generator_test_file.shp not found")
    def test_continuous_map_workflow(self):
        """
        Test the full workflow for continuous map creation.

        Ensures that the continuous map is correctly created with the specified
        parameters and can be rendered without errors.
        """
        # Define input and output file paths
        input_file = os.path.join(
            "tests", "tests_files", "inputs", "html_generator_test_file.shp")
        output_file = os.path.join("tests", "tests_files", "outputs",
                                   "html_generator_InteractiveContinuousMap_test_validation_file.html")
        # Load the test file
        gdf = gpd.read_file(input_file)
        # Initialize map creator with appropriate columns
        map_creator = InteractiveContinuousHtmlMap(
            gdf=gdf,
            color_column='POPULATION',
            tooltip_columns=['ID', 'DEPTO_ID', 'MUN_ID',
                             'MUN_NAME', 'POPULATION', 'REGION']
        )
        # Create the map without saving it
        with mock.patch('folium.Map.save'):
            map_obj = map_creator.create(
                fill_opacity=0.6,
                title='Medellín - Population',
                output_file=output_file
            )
        # Verify the map was created successfully
        assert isinstance(map_obj, folium.Map)
        # Verify it has the expected components
        assert any(isinstance(child, folium.features.GeoJson)
                   for child in map_obj._children.values())
        # Verify there's a colormap for continuous data
        assert any(isinstance(child, cm.LinearColormap)
                   for child in map_obj._children.values())
        # Verify tooltip columns are present
        geo_json_layers = [child for child in map_obj._children.values()
                           if isinstance(child, folium.features.GeoJson)]

    @pytest.mark.skipif(not os.path.exists(os.path.join("tests", "tests_files", "inputs", "html_generator_test_file.shp")),
                        reason="Test file html_generator_test_file.shp not found")
    def test_integrated_map_comparison(self):
        """
        Test comparison between categorical and continuous maps.

        Ensures that both map types can be created from the same dataset
        and verify their differences in visualization approach.
        """
        # Define input file path
        input_file = os.path.join(
            "tests", "tests_files", "inputs", "html_generator_test_file.shp")
        # Load the test file
        gdf = gpd.read_file(input_file)
        # Create both types of maps
        with mock.patch('folium.Map.save'):
            categorical_creator = InteractiveCategoricalHtmlMap(
                gdf=gdf,
                color_column='REGION',
                tooltip_columns=['ID', 'REGION', 'POPULATION']
            )
            categorical_map = categorical_creator.create(title='Categorical Map',
                                                         output_file='html_test_temp_dir/categorical_map.html')
            continuous_creator = InteractiveContinuousHtmlMap(
                gdf=gdf,
                color_column='POPULATION',
                tooltip_columns=['ID', 'REGION', 'POPULATION']
            )
            continuous_map = continuous_creator.create(title='Continuous Map',
                                                       output_file='html_test_temp_dir/continuous_map.html')
            os.rmdir('html_test_temp_dir')
        # Verify the categorical map has the legend HTML injected
        categorical_html = str(categorical_map.get_root())
        # Verify the continuous map has a colormap added
        assert any(isinstance(child, cm.LinearColormap)
                   for child in continuous_map._children.values())
        # Verify different styling mechanisms
        categorical_geo_json = next((child for child in categorical_map._children.values()
                                     if isinstance(child, folium.features.GeoJson)), None)
        continuous_geo_json = next((child for child in continuous_map._children.values()
                                    if isinstance(child, folium.features.GeoJson)), None)
        if categorical_geo_json and continuous_geo_json:
            # Categorical uses discrete colors, continuous uses a gradient
            assert categorical_geo_json.style_function != continuous_geo_json.style_function
