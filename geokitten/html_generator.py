"""
This module provides classes and methods for creating interactive HTML maps using Folium.

It supports both categorical and continuous choropleth maps, with customizable color schemes, 
tooltips, legends, and map styling. The module includes utility classes for validating input 
attributes, generating color palettes, and applying map styles.

Author: Sergio A. Pelayo Escalera
Email: sergioapelayoe@gmail.com
Created: 2025-03-20
Last-modified: 2025-03-24
Version: 0.1.1
"""

__version__ = "0.1.1"
__author__ = "Sergio A. Pelayo Escalera (sergioapelayoe@gmail.com)"
__collabs__ = [""]
__created__ = "2025-03-20"
__last_updated__ = "2025-03-24"

import os
import geopandas as gpd
import pandas as pd
import folium
from folium.features import GeoJsonTooltip
import branca.colormap as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random
from typing import List, Dict, Optional, Tuple
from geokitten.gdf_standardization import StandardGeodataframe


class _MapAttributeValidator:
    """
    Class responsible for validating and processing input attributes for map creation.

    Attributes:
    -----------
    gdf : gpd.GeoDataFrame
        The GeoDataFrame containing polygon geometries and attribute data.
    color_column : str
        The column name used for coloring the polygons.
    tooltip_columns : List[str]
        List of column names to display in the tooltip on hover.
    """

    def __init__(self,
                 gdf: gpd.GeoDataFrame,
                 color_column: str,
                 tooltip_columns: List[str]):
        """
        Initialize the validator with the base required parameters.

        Parameters:
        -----------
        gdf : GeoDataFrame
            The GeoDataFrame containing the polygon geometries and attribute data.
        color_column : str
            The column name to use for coloring the polygons.
        tooltip_columns : list
            List of column names to display in the tooltip on hover.
        """
        self.gdf = self._validate_geodataframe(gdf)
        self.color_column = self._validate_color_column(color_column)
        self.tooltip_columns = self._validate_tooltip_columns(tooltip_columns)

    def _validate_geodataframe(self,
                               gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """
        Validate that the input is a GeoDataFrame.

        Parameters:
        -----------
        gdf : GeoDataFrame
            The input GeoDataFrame.

        Returns:
        --------
        GeoDataFrame
            The validated GeoDataFrame.
        """
        if not isinstance(gdf, gpd.GeoDataFrame):
            raise TypeError("Input must be a GeoDataFrame")
        return gdf

    def _validate_color_column(self,
                               color_column: str) -> str:
        """
        Validate that the color column exists in the GeoDataFrame.

        Parameters:
        -----------
        color_column : str
            The column name to validate.

        Returns:
        --------
        str
            The validated column name.
        """
        if color_column not in self.gdf.columns:
            raise ValueError(
                f"Column '{color_column}' not found in GeoDataFrame")
        return color_column

    def _validate_tooltip_columns(self,
                                  tooltip_columns: List[str]) -> List[str]:
        """
        Validate and filter tooltip columns.

        Parameters:
        -----------
        tooltip_columns : list
            List of column names to validate.

        Returns:
        --------
        list
            List of valid tooltip columns.
        """
        valid_tooltip_columns = [
            col for col in tooltip_columns if col in self.gdf.columns]
        if len(valid_tooltip_columns) == 0:
            raise ValueError(
                "None of the specified tooltip columns exist in the GeoDataFrame")
        if len(valid_tooltip_columns) > 8:
            print(
                f"Warning: Only the first 8 tooltip columns will be used: {valid_tooltip_columns[:8]}")
            valid_tooltip_columns = valid_tooltip_columns[:8]
        return valid_tooltip_columns

    def _get_map_center(self,
                        location: Optional[List[float]] = None) -> List[float]:
        """
        Determine the center location for the map.

        Parameters:
        -----------
        location : list, optional
            Center coordinates [lat, lon] for the map. If None, uses the centroid of the GeoDataFrame.

        Returns:
        --------
        list
            [latitude, longitude] coordinates for map center.
        """
        if location is not None:
            return location
        if self.gdf.crs and not self.gdf.crs.is_geographic:
            centroid = self.gdf.to_crs(epsg=4326).geometry.unary_union.centroid
        else:
            centroid = self.gdf.geometry.unary_union.centroid
        return [centroid.y, centroid.x]

    def _prepare_geodata(self) -> gpd.GeoDataFrame:
        """
        Prepare the geodata for mapping, ensuring it's in geographic coordinates.

        Returns:
        --------
        GeoDataFrame
            The geodata in the correct projection for mapping.
        """
        if self.gdf.crs and not self.gdf.crs.is_geographic:
            return self.gdf.to_crs(epsg=4326)
        return self.gdf

    def _get_bounds(self) -> List[List[float]]:
        """
        Get the bounds of the GeoDataFrame in geographic coordinates.

        Returns:
        --------
        list of lists
            [[min_lat, min_lon], [max_lat, max_lon]] for Folium.
        """
        if self.gdf.crs and not self.gdf.crs.is_geographic:
            bounds = self.gdf.to_crs(epsg=4326).total_bounds
        else:
            bounds = self.gdf.total_bounds
        return [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]


class _CategoricalColorPaletteGenerator:
    """
    Class responsible for generating color palettes for categorical map visualization.

    Attributes:
    -----------
    color_scheme : str
        The color scheme to use (e.g., 'tab20', 'Set1').
    custom_colors : Optional[List[str]]
        List of custom hex color codes, if provided.
    """

    def __init__(self, color_scheme: str = 'tab20'):
        """
        Initialize the categorical color palette generator.

        Parameters:
        -----------
        color_scheme : str
            Color scheme to use (matplotlib colormaps) or 'custom' for custom colors.
        """
        self.color_scheme = color_scheme
        self.custom_colors = None

    def _validate_hex_colors(self,
                             colors: List[str]) -> List[str]:
        """
        Validate that all colors are valid hex codes.

        Parameters:
        -----------
        colors : list
            List of hex color codes to validate.

        Returns:
        --------
        list
            List of valid hex color codes.

        Raises:
        -------
        ValueError: If any color is not a valid hex code.
        """
        for color in colors:
            if not color.startswith('#') or not all(c in '0123456789ABCDEFabcdef' for c in color[1:]):
                raise ValueError(f"Invalid hex color code: {color}")
        return colors

    def _set_custom_colors(self,
                           colors: List[str]) -> None:
        """
        Set custom colors for the palette.

        Parameters:
        -----------
        colors : list
            List of hex color codes to use.

        Raises:
        -------
        ValueError: If any color is not a valid hex code.
        """
        if colors is None:
            return
        colors = self._validate_hex_colors(colors)
        self.custom_colors = colors
        self.color_scheme = 'custom'

    def _generate_distinct_colors(self,
                                  n: int) -> List[str]:
        """
        Generate n visually distinct colors.

        Parameters:
        -----------
        n : int
            Number of distinct colors to generate.

        Returns:
        --------
        list
            List of hex color codes.
        """
        if self.color_scheme == 'custom' and self.custom_colors:
            return self._use_custom_colors(n)
        if n <= 20:
            colors = self._try_specific_colormap(n)
            if colors:
                return colors
        all_colors = self._collect_colors_from_multiple_sources(n)
        if len(all_colors) < n:
            all_colors = self._add_random_colors(all_colors, n)
        return all_colors[:n]

    def _use_custom_colors(self,
                           n: int) -> List[str]:
        """
        Use custom colors with cycling if needed.

        Parameters:
        -----------
        n : int
            Number of colors required.

        Returns:
        --------
        list
            List of hex color codes.
        """
        if len(self.custom_colors) < n:
            return [self.custom_colors[i % len(self.custom_colors)] for i in range(n)]
        return self.custom_colors[:n]

    def _try_specific_colormap(self,
                               n: int) -> Optional[List[str]]:
        """
        Try to use the specific colormap.

        Parameters:
        -----------
        n : int
            Number of colors required.

        Returns:
        --------
        list or None
            List of hex color codes if the colormap is available, otherwise None.
        """
        try:
            colormap = plt.colormaps[self.color_scheme]
            return [mcolors.rgb2hex(colormap(i/(n-1) if n > 1 else 0)) for i in range(n)]
        except (KeyError, ValueError):
            return None

    def _collect_colors_from_multiple_sources(self,
                                              n: int) -> List[str]:
        """
        Collect colors from multiple colormaps and standard color sets.

        Parameters:
        -----------
        n : int
            Number of colors required.

        Returns:
        --------
        list
            List of hex color codes.
        """
        all_colors = []
        self._add_colors_from_colormaps(all_colors, n)
        all_colors.extend(list(mcolors.TABLEAU_COLORS.values()))
        if n > len(all_colors):
            all_colors.extend(list(mcolors.CSS4_COLORS.values()))
        return list(dict.fromkeys(all_colors))

    def _add_colors_from_colormaps(self,
                                   all_colors: List[str], n: int) -> None:
        """
        Add colors from various colormaps.

        Parameters:
        -----------
        all_colors : list
            List to append colors to.
        n : int
            Number of colors required.
        """
        for cmap_name in ['Set1', 'Pastel1', 'Paired', 'Dark2', 'tab10', 'tab20']:
            try:
                cmap = plt.colormaps[cmap_name]
                steps = min(20, n)
                all_colors.extend(
                    [mcolors.rgb2hex(cmap(i/(steps-1) if steps > 1 else 0)) for i in range(steps)])
            except (KeyError, ValueError):
                pass

    def _add_random_colors(self,
                           existing_colors: List[str], n: int) -> List[str]:
        """
        Generate random colors with good contrast against existing ones.

        Parameters:
        -----------
        existing_colors : list
            List of existing colors.
        n : int
            Number of colors required.

        Returns:
        --------
        list
            List of hex color codes.
        """
        all_colors = existing_colors.copy()
        existing_hsv = [mcolors.rgb_to_hsv(
            mcolors.to_rgb(color)) for color in all_colors]
        while len(all_colors) < n:
            new_color = self._generate_distinct_random_color(existing_hsv)
            if new_color:
                all_colors.append(new_color)
        return all_colors

    def _generate_distinct_random_color(self,
                                        existing_hsv: List[List[float]]) -> Optional[str]:
        """
        Generate a single random color that is distinct from existing colors.

        Parameters:
        -----------
        existing_hsv : list
            List of HSV values of existing colors.

        Returns:
        --------
        str or None
            A distinct hex color code, or None if no distinct color could be generated.
        """
        h = random.random()
        s = 0.5 + random.random() * 0.5  # 0.5-1.0 saturation
        v = 0.9 + random.random() * 0.1  # 0.9-1.0 value/brightness
        if self._is_color_distinct(h, s, v, existing_hsv):
            r, g, b = mcolors.hsv_to_rgb([h, s, v])
            new_color = mcolors.rgb2hex((r, g, b))
            existing_hsv.append([h, s, v])
            return new_color
        return None

    def _is_color_distinct(self,
                           h: float,
                           s: float,
                           v: float,
                           existing_hsv: List[List[float]]) -> bool:
        """
        Check if a color is distinct from all existing colors.
        Compute the distance in HSV space (with wrap-around for hue).

        Parameters:
        -----------
        h : float
            Hue value of the color.
        s : float
            Saturation value of the color.
        v : float
            Value (brightness) of the color.
        existing_hsv : list
            List of HSV values of existing colors.

        Returns:
        --------
        bool
            True if the color is sufficiently distinct, False otherwise.
        """
        min_distance = 1.0
        for existing in existing_hsv:
            h_diff = min(abs(h - existing[0]), 1 - abs(h - existing[0]))
            s_diff = abs(s - existing[1])
            v_diff = abs(v - existing[2])
            distance = (h_diff**2 + s_diff**2 + v_diff**2)**0.5
            min_distance = min(min_distance, distance)
        return min_distance > 0.15

    def _create_categorical_color_dict(self,
                                       values: List) -> Dict:
        """
        Create a dictionary mapping values to colors.

        Parameters:
        -----------
        values : list
            List of unique values to map to colors.

        Returns:
        --------
        dict
            Dictionary mapping values to hex color codes.
        """
        unique_values = sorted(values) if all(isinstance(
            x, (int, float)) for x in values) else values
        if self.color_scheme == 'custom' and self.custom_colors and len(self.custom_colors) == 1:
            return {value: self.custom_colors[0] for value in unique_values}
        colors = self._generate_distinct_colors(len(unique_values))
        return {value: colors[i] for i, value in enumerate(unique_values)}


class _ContinuousColorPaletteGenerator:
    """
    Class responsible for generating color palettes for continuous map visualization.

    Attributes:
    -----------
    color_scheme : str
        The color scheme to use (e.g., 'viridis', 'coolwarm').
    """

    def __init__(self, color_scheme: str = 'viridis'):
        """
        Initialize the continuous color palette generator.

        Parameters:
        -----------
        color_scheme : str
            Color scheme to use from matplotlib colormaps.
        """
        self.color_scheme = color_scheme
        self._validate_color_scheme()

    def _validate_color_scheme(self) -> None:
        """
        Validate that the specified color scheme exists in Matplotlib.

        This method checks if the `color_scheme` attribute is a valid Matplotlib colormap. 
        If the color scheme is not found, it raises a `ValueError` and provides a list 
        of available colormaps.

        Raises:
        -------
        ValueError:
            If the specified color scheme is not found in Matplotlib colormaps.
        """
        if self.color_scheme not in plt.colormaps():
            available_cmaps = ", ".join(sorted([cm for cm in plt.colormaps()
                                                if not cm.endswith('_r') and not cm.startswith('_')]))
            raise ValueError(f"Color scheme '{self.color_scheme}' not found in matplotlib. "
                             f"Available color schemes include: {available_cmaps}")

    def _create_continuous_colormap(self,
                                    min_value: float,
                                    max_value: float) -> cm.LinearColormap:
        """
        Create a continuous colormap for numeric data.
        Get colors from matplotlib colormap

        Parameters:
        -----------
        min_value : float
            Minimum value in the data
        max_value : float
            Maximum value in the data

        Returns:
        --------
        branca.colormap.LinearColormap
            A colormap object for the data range
        """
        cmap = plt.get_cmap(self.color_scheme)
        colors = [mcolors.rgb2hex(cmap(i/255)) for i in range(256)]
        return cm.LinearColormap(
            colors=colors,
            vmin=min_value,
            vmax=max_value
        )


class _MapStyleGenerator:
    """
    Class responsible for generating the style for the map.

    Attributes:
    -----------
    fill_opacity : float
        Opacity of the polygon fill (0-1).
    line_opacity : float
        Opacity of the polygon borders (0-1).
    line_weight : int
        Width of the polygon borders.
    """

    def __init__(self,
                 fill_opacity: float = 0.5,
                 line_opacity: float = 1.0,
                 line_weight: int = 1):
        """
        Initialize the style generator.

        Parameters:
        -----------
        fill_opacity : float
            Opacity of the polygon fill (0-1)
        line_opacity : float
            Opacity of the polygon borders (0-1)
        line_weight : int
            Width of the polygon borders
        """
        self.fill_opacity = fill_opacity
        self.line_opacity = line_opacity
        self.line_weight = line_weight

    def _create_categorical_style(self,
                                  color_column: str,
                                  color_dict: Dict) -> Tuple[callable, str]:
        """
        Create style function and legend HTML for categorical data.

        Parameters:
        -----------
        color_column : str
            The column name used for coloring
        color_dict : dict
            Dictionary mapping values to colors

        Returns:
        --------
        tuple
            (style_function, legend_html)
        """
        def style_function(x): return {
            'fillColor': color_dict.get(x['properties'][color_column], '#FFFFFF'),
            'color': 'black',
            'weight': self.line_weight,
            'fillOpacity': self.fill_opacity,
            'opacity': self.line_opacity
        }
        legend_html = '''
        <div style="position: fixed; bottom: 50px; left: 50px; z-index: 1000; padding: 10px; 
        background-color: white; border-radius: 5px; border: 2px solid grey; opacity: 0.8; 
        max-height: 300px; overflow-y: auto;">
        <p style="text-align: center; margin-bottom: 5px;"><strong>{}</strong></p>
        '''.format(color_column)
        for value, color in color_dict.items():
            legend_html += '''
            <div style="display: flex; align-items: center; margin: 3px;">
            <div style="width: 15px; height: 15px; background-color: {}; margin-right: 5px;"></div>
            <span>{}</span>
            </div>
            '''.format(color, value)
        legend_html += '</div>'
        return style_function, legend_html

    def _create_continuous_style(self,
                                 color_column: str,
                                 colormap: cm.LinearColormap) -> callable:
        """
        Create style function for continuous data.

        Parameters:
        -----------
        color_column : str
            The column name used for coloring
        colormap : branca.colormap.LinearColormap
            The colormap to use

        Returns:
        --------
        callable
            Style function for the GeoJSON
        """
        return lambda x: {
            'fillColor': colormap(x['properties'][color_column]),
            'color': 'black',
            'weight': self.line_weight,
            'fillOpacity': self.fill_opacity,
            'opacity': self.line_opacity
        }


class _InteractiveMapCreator:
    """
    Class responsible for creating the interactive map.

    Attributes:
    -----------
    validator : _MapAttributeValidator
        Validator object with validated parameters.
    """

    def __init__(self,
                 validator: _MapAttributeValidator):
        """
        Initialize the map creator.

        Parameters:
        -----------
        validator : _MapAttributeValidator
            Validator object with validated parameters
        """
        self.validator = validator

    def _create_map(self,
                    zoom_start: Optional[int] = None,
                    location: Optional[List[float]] = None) -> folium.Map:
        """
        Create the base folium map.

        Parameters:
        -----------
        zoom_start : int, optional
            Initial zoom level (if None, will be auto-calculated based on bounds)
        location : list, optional
            Center coordinates [lat, lon]

        Returns:
        --------
        folium.Map
            Base map object
        """
        center = self.validator._get_map_center(location)
        if zoom_start is None:
            m = folium.Map(location=center)
            bounds = self.validator._get_bounds()
            m.fit_bounds(bounds)
        else:
            m = folium.Map(location=center, zoom_start=zoom_start)
        return m

    def _add_choropleth_layer(self,
                              m: folium.Map,
                              style_function: callable,
                              tooltip_columns: List[str]) -> folium.Map:
        """
        Add the choropleth layer to the map.

        Parameters:
        -----------
        m : folium.Map
            Base map object
        style_function : callable
            Function that defines the style for each feature
        tooltip_columns : list
            List of column names to display in tooltips

        Returns:
        --------
        folium.Map
            Map with the added layer
        """
        geo_json_data = self.validator._prepare_geodata()
        tooltip = GeoJsonTooltip(
            fields=tooltip_columns,
            aliases=tooltip_columns,
            style=("background-color: white; color: #333333; font-family: arial; font-size: 12px; "
                   "padding: 10px; border-radius: 3px; box-shadow: 3px 3px 3px rgba(0,0,0,0.2);")
        )
        folium.GeoJson(
            data=geo_json_data,
            style_function=style_function,
            tooltip=tooltip,
            name='Polygons'
        ).add_to(m)
        return m

    def _add_controls(self,
                      m: folium.Map) -> folium.Map:
        """
        Add control elements to the map.

        Parameters:
        -----------
        m : folium.Map
            The map object

        Returns:
        --------
        folium.Map
            Map with added controls
        """
        folium.LayerControl().add_to(m)
        return m

    def _add_title(self,
                   m: folium.Map,
                   title: str) -> folium.Map:
        """
        Add a title to the map.

        Parameters:
        -----------
        m : folium.Map
            The map object
        title : str
            The title text

        Returns:
        --------
        folium.Map
            Map with added title
        """
        if title:
            title_html = '''
            <div style="position: fixed; top: 10px; left: 50%; transform: translateX(-50%); z-index: 1000; 
            background-color: white; padding: 10px; border-radius: 5px; border: 2px solid grey; opacity: 0.8;">
            <h3 style="margin: 0;">{}</h3>
            </div>
            '''.format(title)
            m.get_root().html.add_child(folium.Element(title_html))
        return m

    def _add_categorical_legend(self,
                                m: folium.Map,
                                legend_html: str) -> folium.Map:
        """
        Add a categorical legend to the map.

        Parameters:
        -----------
        m : folium.Map
            The map object
        legend_html : str
            HTML for the legend

        Returns:
        --------
        folium.Map
            Map with added legend
        """
        m.get_root().html.add_child(folium.Element(legend_html))
        return m

    def _add_continuous_colormap(self,
                                 m: folium.Map,
                                 colormap: cm.LinearColormap,
                                 legend_name: str) -> folium.Map:
        """
        Add a colormap to the map for continuous data with an improved layout.
        Add the colormap directly to the map first

        Parameters:
        -----------
        m : folium.Map
            The map object
        colormap : branca.colormap.LinearColormap
            The colormap to add
        legend_name : str
            Name for the legend

        Returns:
        --------
        folium.Map
            Map with added colormap
        """
        colormap.caption = legend_name
        colormap.add_to(m)
        return m


class InteractiveCategoricalHtmlMap:
    """
    Class for creating categorical interactive choropleth maps.

    Attributes:
    -----------
    validator : _MapAttributeValidator
        Validator object for input attributes.
    map_creator : _InteractiveMapCreator
        Map creator object for managing the map.
    """

    def __init__(self,
                 gdf: gpd.GeoDataFrame,
                 color_column: str,
                 tooltip_columns: List[str]):
        """
        Create a categorical interactive map.

        This method generates a categorical choropleth map using the provided GeoDataFrame, 
        color column, and styling options. The map is saved as an HTML file.

        Parameters:
        -----------
        color_scheme : str
            Color scheme to use (e.g., 'tab20', 'Set1').
        set_custom_colors : list, optional
            List of hex color codes to use. If only one color is provided, all geometries 
            will be colored the same.
        fill_opacity : float
            Opacity of the polygon fill (0-1).
        line_opacity : float
            Opacity of the polygon borders (0-1).
        line_weight : int
            Width of the polygon borders.
        location : list, optional
            Center coordinates [lat, lon] for the map. If None, the center is auto-calculated.
        zoom_start : int, optional
            Initial zoom level for the map. If None, the zoom level is auto-calculated.
        title : str, optional
            Title for the map. If None, no title is added.
        legend_name : str, optional
            Name for the legend. If None, the column name used for coloring is used as the legend name.
        output_file : str
            Filename to save the HTML map.

        Returns:
        --------
        folium.Map
            The created map object.

        Examples:
        ---------
        Example 1: Creating a map with default settings
        -----------------------------------------------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Polygon
        >>> from src.html_generator import InteractiveCategoricalHtmlMap
        >>> data = {
        ...     'geometry': [Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
        ...     'region': ['Region A']
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
        >>> map_creator = InteractiveCategoricalHtmlMap(
        ...     gdf=gdf, 
        ...     color_column='region', 
        ...     tooltip_columns=['region']
        ... )
        >>> map_creator.create(output_file='categorical_map.html')

        Example 2: Customizing the map with a different color scheme and title
        ----------------------------------------------------------------------
        >>> map_creator.create(
        ...     color_scheme='Set1',
        ...     fill_opacity=0.7,
        ...     line_opacity=0.8,
        ...     line_weight=2,
        ...     title='Regional Map',
        ...     legend_name='Regions',
        ...     output_file='custom_categorical_map.html'
        ... )

        Example 3: Using custom colors for the map
        ------------------------------------------
        >>> map_creator.create(
        ...     set_custom_colors=['#FF0000', '#00FF00', '#0000FF'],
        ...     title='Custom Color Map',
        ...     output_file='custom_color_map.html'
        ... )

        Example 4: Specifying a custom map center and zoom level
        --------------------------------------------------------
        >>> map_creator.create(
        ...     location=[0.5, 0.5],
        ...     zoom_start=10,
        ...     output_file='centered_categorical_map.html'
        ... )
        """
        self.validator = _MapAttributeValidator(
            StandardGeodataframe(gdf, crs=gdf.crs, remove_geni=False),
            color_column,
            tooltip_columns)
        self.map_creator = _InteractiveMapCreator(self.validator)

    def _initialize_generators(self,
                               color_scheme: str = 'tab20',
                               set_custom_colors: Optional[List[str]] = None,
                               fill_opacity: float = 0.5,
                               line_opacity: float = 1.0,
                               line_weight: int = 1) -> Tuple[_CategoricalColorPaletteGenerator, _MapStyleGenerator]:
        """
        Initialize the color and style generators.

        Parameters:
        -----------
        color_scheme : str
            Color scheme to use for the map (e.g., 'tab20', 'Set1').
        set_custom_colors : list, optional
            List of custom hex color codes to override the default color scheme.
        fill_opacity : float
            Opacity of the polygon fill (range: 0-1).
        line_opacity : float
            Opacity of the polygon borders (range: 0-1).
        line_weight : int
            Width of the polygon borders.

        Returns:
        --------
        Tuple[_CategoricalColorPaletteGenerator, _MapStyleGenerator]
            The initialized color and style generators.
        """
        color_generator = _CategoricalColorPaletteGenerator(color_scheme)
        if set_custom_colors:
            color_generator._set_custom_colors(set_custom_colors)
        style_generator = _MapStyleGenerator(
            fill_opacity=fill_opacity,
            line_opacity=line_opacity,
            line_weight=line_weight
        )
        return color_generator, style_generator

    def _create_map(self,
                    location: Optional[List[float]] = None,
                    zoom_start: Optional[int] = None) -> folium.Map:
        """
        Create the base folium map.

        Parameters:
        -----------
        location : list, optional
            Center coordinates [lat, lon] for the map. If None, the center is auto-calculated.
        zoom_start : int, optional
            Initial zoom level for the map. If None, the zoom level is auto-calculated.

        Returns:
        --------
        folium.Map
            The base map object.
        """
        return self.map_creator._create_map(zoom_start, location)

    def _apply_categorical_styling(self,
                                   color_generator: _CategoricalColorPaletteGenerator,
                                   style_generator: _MapStyleGenerator,
                                   legend_name: Optional[str] = None) -> Tuple[callable, str]:
        """
        Apply categorical styling to the map.

        Parameters:
        -----------
        color_generator : _CategoricalColorPaletteGenerator
            The color generator used to create a color palette for categorical data.
        style_generator : _MapStyleGenerator
            The style generator used to define the map's visual style.
        legend_name : str, optional
            The name of the legend. If None, the column name used for coloring is used as the legend name.

        Returns:
        --------
        Tuple[callable, str]
            A tuple containing:
            - `style_function`: A callable function that defines the style for each feature.
            - `legend_html`: A string containing the HTML for the legend.
        """
        gdf = self.validator.gdf
        color_column = self.validator.color_column
        legend_name = legend_name or color_column
        unique_values = gdf[color_column].unique()
        color_dict = color_generator._create_categorical_color_dict(
            unique_values)
        style_function, legend_html = style_generator._create_categorical_style(
            color_column, color_dict
        )
        return style_function, legend_html

    def _add_layers_controls(self,
                             m: folium.Map,
                             style_function: callable,
                             legend_html: str,
                             title: Optional[str] = None) -> folium.Map:
        """
        Add layers, controls, and title to the map.

        Parameters:
        -----------
        m : folium.Map
            The base map object.
        style_function : callable
            The style function for the map.
        legend_html : str
            The HTML string for the legend.
        title : str, optional
            The title of the map. If None, no title is added.

        Returns:
        --------
        folium.Map
            The updated map object with added layers, controls, and title.
        """
        gdf = self.validator.gdf
        tooltip_columns = self.validator.tooltip_columns
        m = self.map_creator._add_choropleth_layer(
            m, style_function, tooltip_columns)
        m = self.map_creator._add_categorical_legend(m, legend_html)
        m = self.map_creator._add_controls(m)
        m = self.map_creator._add_title(m, title)
        return m

    def create(self,
               color_scheme: str = 'tab20',
               set_custom_colors: Optional[List[str]] = None,
               fill_opacity: float = 0.5,
               line_opacity: float = 1.0,
               line_weight: int = 1,
               location: Optional[List[float]] = None,
               zoom_start: Optional[int] = None,
               title: Optional[str] = None,
               output_file: str = 'categorical_map.html') -> folium.Map:
        """
        Create a categorical interactive map.

        This method generates a categorical choropleth map using the provided GeoDataFrame, 
        color column, and styling options. The map is saved as an HTML file.

        Parameters:
        -----------
        color_scheme : str
            Color scheme to use (e.g., 'tab20', 'Set1').
        set_custom_colors : list, optional
            List of hex color codes to use. If only one color is provided, all geometries 
            will be colored the same.
        fill_opacity : float
            Opacity of the polygon fill (0-1).
        line_opacity : float
            Opacity of the polygon borders (0-1).
        line_weight : int
            Width of the polygon borders.
        location : list, optional
            Center coordinates [lat, lon] for the map. If None, the center is auto-calculated.
        zoom_start : int, optional
            Initial zoom level for the map. If None, the zoom level is auto-calculated.
        title : str, optional
            Title for the map. If None, no title is added.
        legend_name : str, optional
            Name for the legend. If None, the column name used for coloring is used as the legend name.
        output_file : str
            Filename to save the HTML map.

        Returns:
        --------
        folium.Map
            The created map object.

        Examples:
        ---------
        Example 1: Creating a map with default settings
        -----------------------------------------------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Polygon
        >>> from src.html_generator import InteractiveCategoricalHtmlMap
        >>> data = {
        ...     'geometry': [Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
        ...     'region': ['Region A']
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
        >>> map_creator = InteractiveCategoricalHtmlMap(
        ...     gdf=gdf, 
        ...     color_column='region', 
        ...     tooltip_columns=['region']
        ... )
        >>> map_creator.create(output_file='categorical_map.html')

        Example 2: Customizing the map with a different color scheme and title
        ----------------------------------------------------------------------
        >>> map_creator.create(
        ...     color_scheme='Set1',
        ...     fill_opacity=0.7,
        ...     line_opacity=0.8,
        ...     line_weight=2,
        ...     title='Regional Map',
        ...     legend_name='Regions',
        ...     output_file='custom_categorical_map.html'
        ... )

        Example 3: Using custom colors for the map
        ------------------------------------------
        >>> map_creator.create(
        ...     set_custom_colors=['#FF0000', '#00FF00', '#0000FF'],
        ...     title='Custom Color Map',
        ...     output_file='custom_color_map.html'
        ... )

        Example 4: Specifying a custom map center and zoom level
        --------------------------------------------------------
        >>> map_creator.create(
        ...     location=[0.5, 0.5],
        ...     zoom_start=10,
        ...     output_file='centered_categorical_map.html'
        ... )
        """
        legend_name = None
        color_generator, style_generator = self._initialize_generators(
            color_scheme, set_custom_colors, fill_opacity, line_opacity, line_weight
        )
        m = self._create_map(location, zoom_start)
        style_function, legend_html = self._apply_categorical_styling(
            color_generator, style_generator, legend_name
        )
        m = self._add_layers_controls(m, style_function, legend_html, title)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        m.save(output_file)
        return m


class InteractiveContinuousHtmlMap:
    """
    Class for creating continuous interactive choropleth maps.

    Attributes:
    -----------
    validator : _MapAttributeValidator
        Validator object for input attributes.
    map_creator : _InteractiveMapCreator
        Map creator object for managing the map.
    """

    def __init__(self,
                 gdf: gpd.GeoDataFrame,
                 color_column: str,
                 tooltip_columns: List[str]):
        """
        Initialize the continuous interactive map creator.

        This constructor sets up the instance by validating the input GeoDataFrame, 
        the color column, and the tooltip columns. It also ensures that the color 
        column is numeric, as required for continuous maps.

        Parameters:
        -----------
        gdf : GeoDataFrame
            The GeoDataFrame containing the polygon geometries and attribute data.
        color_column : str
            The column name to use for coloring the polygons (must be numeric).
        tooltip_columns : list
            List of column names to display in the tooltip on hover.

        Raises:
        -------
        ValueError:
            If the specified color column is not numeric.

        Examples:
        ---------
        Example 1: Initializing with a valid GeoDataFrame
        -------------------------------------------------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Polygon
        >>> from src.html_generator import InteractiveContinuousHtmlMap
        >>> data = {
        ...     'geometry': [Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
        ...     'population_density': [1000],
        ...     'region_name': ['Region A']
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
        >>> map_creator = InteractiveContinuousHtmlMap(
        ...     gdf=gdf, 
        ...     color_column='population_density', 
        ...     tooltip_columns=['region_name', 'population_density']
        ... )

        Example 2: Invalid color column (non-numeric)
        ---------------------------------------------
        >>> data = {
        ...     'geometry': [Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
        ...     'region_name': ['Region A']
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
        >>> map_creator = InteractiveContinuousHtmlMap(
        ...     gdf=gdf, 
        ...     color_column='region_name', 
        ...     tooltip_columns=['region_name']
        ... )
        ValueError: Color column 'region_name' must be numeric for continuous maps
        """
        self.validator = _MapAttributeValidator(
            StandardGeodataframe(gdf, crs=gdf.crs, remove_geni=False),
            color_column,
            tooltip_columns)
        self.map_creator = _InteractiveMapCreator(self.validator)
        self._validate_color_column(gdf, color_column)

    def _validate_color_column(self,
                               gdf: gpd.GeoDataFrame,
                               color_column: str) -> None:
        """
        Validate that the specified color column is numeric.

        This method checks if the `color_column` in the provided GeoDataFrame is of a numeric data type.
        If the column is not numeric, it raises a `ValueError`.

        Parameters:
        -----------
        gdf : GeoDataFrame
            The GeoDataFrame containing the data.
        color_column : str
            The name of the column to validate.

        Raises:
        -------
        ValueError:
            If the specified column is not numeric.
        """
        if not pd.api.types.is_numeric_dtype(gdf[color_column]):
            raise ValueError(
                f"Color column '{color_column}' must be numeric for continuous maps")

    def _initialize_generators(self,
                               color_scheme: str,
                               fill_opacity: float,
                               line_opacity: float,
                               line_weight: int) -> Tuple[_ContinuousColorPaletteGenerator, _MapStyleGenerator]:
        """
        Initialize the color and style generators.

        Parameters:
        -----------
        color_scheme : str
            Matplotlib colormap name to use
        fill_opacity : float
            Opacity of the polygon fill (0-1)
        line_opacity : float
            Opacity of the polygon borders (0-1)
        line_weight : int
            Width of the polygon borders

        Returns:
        --------
        Tuple[_ContinuousColorPaletteGenerator, _MapStyleGenerator]
            The initialized color and style generators
        """
        color_generator = _ContinuousColorPaletteGenerator(color_scheme)
        style_generator = _MapStyleGenerator(
            fill_opacity=fill_opacity,
            line_opacity=line_opacity,
            line_weight=line_weight
        )
        return color_generator, style_generator

    def _create_map(self,
                    location: Optional[List[float]],
                    zoom_start: Optional[int]) -> folium.Map:
        """
        Create the base folium map.

        Parameters:
        -----------
        location : list, optional
            Center coordinates [lat, lon] for the map
        zoom_start : int, optional
            Initial zoom level (if None, will auto-fit to data)

        Returns:
        --------
        folium.Map
            The base map object
        """
        return self.map_creator._create_map(zoom_start, location)

    def _apply_continuous_styling(self,
                                  color_generator: _ContinuousColorPaletteGenerator,
                                  style_generator: _MapStyleGenerator,
                                  legend_name: Optional[str]) -> Tuple[callable, cm.LinearColormap]:
        """
        Apply continuous styling to the map.

        Parameters:
        -----------
        color_generator : _ContinuousColorPaletteGenerator
            The color generator for continuous data
        style_generator : _MapStyleGenerator
            The style generator for the map
        legend_name : str, optional
            Name for the legend

        Returns:
        --------
        Tuple[callable, cm.LinearColormap]
            The style function and colormap for the map
        """
        gdf = self.validator.gdf
        color_column = self.validator.color_column
        legend_name = legend_name or color_column
        min_value = gdf[color_column].min()
        max_value = gdf[color_column].max()
        colormap = color_generator._create_continuous_colormap(
            min_value, max_value)
        style_function = style_generator._create_continuous_style(
            color_column, colormap)
        return style_function, colormap

    def _add_layers_controls(self,
                             m: folium.Map,
                             style_function: callable,
                             colormap: cm.LinearColormap,
                             legend_name: str,
                             title: Optional[str]) -> folium.Map:
        """
        Add layers, controls, and title to the map.

        Parameters:
        -----------
        m : folium.Map
            The base map object
        style_function : callable
            The style function for the map
        colormap : cm.LinearColormap
            The colormap for continuous data
        legend_name : str
            Name for the legend
        title : str, optional
            Title for the map

        Returns:
        --------
        folium.Map
            The updated map object with layers and controls
        """
        gdf = self.validator.gdf
        tooltip_columns = self.validator.tooltip_columns
        m = self.map_creator._add_choropleth_layer(
            m, style_function, tooltip_columns)
        m = self.map_creator._add_continuous_colormap(m, colormap, legend_name)
        m = self.map_creator._add_controls(m)
        m = self.map_creator._add_title(m, title)
        return m

    def create(self,
               color_scheme: str = 'viridis',
               fill_opacity: float = 0.5,
               line_opacity: float = 1.0,
               line_weight: int = 1,
               location: Optional[List[float]] = None,
               zoom_start: Optional[int] = None,
               title: Optional[str] = None,
               legend_name: Optional[str] = None,
               output_file: str = 'continuous_map.html') -> folium.Map:
        """
        Create a continuous interactive map.

        This method generates a continuous choropleth map using the provided GeoDataFrame, 
        color column, and styling options. The map is saved as an HTML file.

        Parameters:
        -----------
        color_scheme : str
            Matplotlib colormap name to use (e.g., 'viridis', 'coolwarm', 'Spectral', 'RdBu').
        fill_opacity : float
            Opacity of the polygon fill (0-1).
        line_opacity : float
            Opacity of the polygon borders (0-1).
        line_weight : int
            Width of the polygon borders.
        location : list, optional
            Center coordinates [lat, lon] for the map. If None, the center is auto-calculated.
        zoom_start : int, optional
            Initial zoom level for the map. If None, the zoom level is auto-calculated.
        title : str, optional
            Title for the map. If None, no title is added.
        legend_name : str, optional
            Name for the legend. If None, the column name used for coloring is used as the legend name.
        output_file : str
            Filename to save the HTML map.

        Returns:
        --------
        folium.Map
            The created map object.

        Examples:
        ---------
        Example 1: Creating a map with default settings
        -----------------------------------------------
        >>> import geopandas as gpd
        >>> from shapely.geometry import Polygon
        >>> from src.html_generator import InteractiveContinuousHtmlMap
        >>> data = {
        ...     'geometry': [Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
        ...     'population_density': [1000]
        ... }
        >>> gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
        >>> map_creator = InteractiveContinuousHtmlMap(
        ...     gdf=gdf, 
        ...     color_column='population_density', 
        ...     tooltip_columns=['population_density']
        ... )
        >>> map_creator.create(output_file='population_density_map.html')

        Example 2: Customizing the map with a different colormap and title
        ------------------------------------------------------------------
        >>> map_creator.create(
        ...     color_scheme='coolwarm',
        ...     fill_opacity=0.7,
        ...     line_opacity=0.8,
        ...     line_weight=2,
        ...     title='Population Density Map',
        ...     legend_name='Density (people per km²)',
        ...     output_file='custom_population_density_map.html'
        ... )

        Example 3: Specifying a custom map center and zoom level
        --------------------------------------------------------
        >>> map_creator.create(
        ...     location=[0.5, 0.5],
        ...     zoom_start=10,
        ...     output_file='centered_population_density_map.html'
        ... )
        """
        color_generator, style_generator = self._initialize_generators(
            color_scheme, fill_opacity, line_opacity, line_weight
        )
        m = self._create_map(location, zoom_start)
        style_function, colormap = self._apply_continuous_styling(
            color_generator, style_generator, legend_name
        )
        m = self._add_layers_controls(
            m, style_function, colormap, legend_name, title)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        m.save(output_file)
        return m
