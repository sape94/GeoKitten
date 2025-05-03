# HTML Generator Module

A Python module for creating interactive HTML maps using Folium. The module supports both categorical and continuous choropleth maps with customizable color schemes, tooltips, legends, and map styling.

- **Author**: Sergio A. Pelayo Escalera
- **Email**: sergioapelayoe@gmail.com
- **Version**: 0.1.0
- **Created**: 2025-04-09
- **Last Modified**: 2025-04-09

## Overview

This module provides classes and methods for generating interactive maps in HTML format. It includes utilities for validating input attributes, generating color palettes, and applying map styles.

### Key Features:

- Categorical and continuous choropleth maps.
- Customizable color schemes and styles.
- Tooltip and legend support.
- Integration with GeoDataFrames for spatial data visualization.

## Requirements

The module has the following dependencies:

- `geopandas`
- `folium`
- `branca`
- `matplotlib`
- `shapely`

## Classes

### InteractiveCategoricalHtmlMap

A class for creating categorical interactive choropleth maps.

#### Initialization

```python
from src.html_generator import InteractiveCategoricalHtmlMap

# Initialize with a GeoDataFrame
map_creator = InteractiveCategoricalHtmlMap(
    gdf=gdf,
    color_column="region",
    tooltip_columns=["region", "population"]
)
```

#### Methods

##### `create(color_scheme='tab20', set_custom_colors=None, fill_opacity=0.5, line_opacity=1.0, line_weight=1, location=None, zoom_start=None, title=None, legend_name=None, output_file='categorical_map.html')`

Generates a categorical choropleth map and saves it as an HTML file.

**Parameters**:

- `color_scheme` (str, optional): Color scheme to use (e.g., 'tab20', 'Set1'). Defaults to 'tab20'.
- `set_custom_colors` (list, optional): List of custom hex color codes. Defaults to None.
- `fill_opacity` (float, optional): Opacity of the polygon fill (0-1). Defaults to 0.5.
- `line_opacity` (float, optional): Opacity of the polygon borders (0-1). Defaults to 1.0.
- `line_weight` (int, optional): Width of the polygon borders. Defaults to 1.
- `location` (list, optional): Center coordinates [lat, lon] for the map. Defaults to None.
- `zoom_start` (int, optional): Initial zoom level for the map. Defaults to None.
- `title` (str, optional): Title for the map. Defaults to None.
- `output_file` (str, optional): Filename to save the HTML map. Defaults to 'categorical_map.html'.

**Returns**:

- `folium.Map`: The created map object.

**Example**:

```python
map_creator = InteractiveCategoricalHtmlMap(
    gdf=gdf,
    color_column="region",
    tooltip_columns=["region", "population"]
)
map_creator.create(
    color_scheme="Set1",
    fill_opacity=0.7,
    title="Regional Map",
    output_file="categorical_map.html"
)
```

### InteractiveContinuousHtmlMap

A class for creating continuous interactive choropleth maps.

#### Initialization

```python
from src.html_generator import InteractiveContinuousHtmlMap

# Initialize with a GeoDataFrame
map_creator = InteractiveContinuousHtmlMap(
    gdf=gdf,
    color_column="population_density",
    tooltip_columns=["region", "population_density"]
)
```

#### Methods

##### `create(color_scheme='viridis', fill_opacity=0.5, line_opacity=1.0, line_weight=1, location=None, zoom_start=None, title=None, legend_name=None, output_file='continuous_map.html')`

Generates a continuous choropleth map and saves it as an HTML file.

**Parameters**:

- `color_scheme` (str, optional): Matplotlib colormap name to use (e.g., 'viridis', 'coolwarm'). Defaults to 'viridis'.
- `fill_opacity` (float, optional): Opacity of the polygon fill (0-1). Defaults to 0.5.
- `line_opacity` (float, optional): Opacity of the polygon borders (0-1). Defaults to 1.0.
- `line_weight` (int, optional): Width of the polygon borders. Defaults to 1.
- `location` (list, optional): Center coordinates [lat, lon] for the map. Defaults to None.
- `zoom_start` (int, optional): Initial zoom level for the map. Defaults to None.
- `title` (str, optional): Title for the map. Defaults to None.
- `legend_name` (str, optional): Name for the legend. Defaults to None.
- `output_file` (str, optional): Filename to save the HTML map. Defaults to 'continuous_map.html'.

**Example**:

```python
map_creator = InteractiveContinuousHtmlMap(
    gdf=gdf,
    color_column="population_density",
    tooltip_columns=["region", "population_density"]
)
map_creator.create(
    color_scheme="coolwarm",
    fill_opacity=0.7,
    title="Population Density Map",
    output_file="population_density_map.html"
)
```

## Usage Examples

### Creating a Categorical Map

```python
from src.html_generator import InteractiveCategoricalHtmlMap

# Create a GeoDataFrame
import geopandas as gpd
from shapely.geometry import Polygon

data = {
    "geometry": [Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
    "region": ["Region A"]
}
gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")

# Generate a categorical map
map_creator = InteractiveCategoricalHtmlMap(
    gdf=gdf,
    color_column="region",
    tooltip_columns=["region"]
)
map_creator.create(output_file="categorical_map.html")
```

### Creating a Continuous Map

```python
from src.html_generator import InteractiveContinuousHtmlMap

# Create a GeoDataFrame
import geopandas as gpd
from shapely.geometry import Polygon

data = {
    "geometry": [Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])],
    "population_density": [1000]
}
gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")

# Generate a continuous map
map_creator = InteractiveContinuousHtmlMap(
    gdf=gdf,
    color_column="population_density",
    tooltip_columns=["population_density"]
)
map_creator.create(output_file="population_density_map.html")
```

## Notes

- All GeoDataFrames must use the WGS 84 (EPSG:4326) coordinate system.
- The module supports both Polygon and MultiPolygon geometries.
- Maps are saved as standalone HTML files that can be opened in any web browser.
