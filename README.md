# GeoKitten

A Python package for standardizing, converting, and visualizing geospatial data with an emphasis on simplicity and extensibility.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
- [Module Overview](#module-overview)
  - [1. StandardGeodataframe Class](#1-standardgeodataframe-class)
    - [Constructor](#constructor)
    - [Class Methods](#class-methods)
    - [Methods](#methods)
  - [2. KML Converter Classes](#2-kml-converter-classes)
    - [KMLsToGeodataframe](#kmlstogeodataframe)
    - [GeodataframeToKMLs](#geodataframetokml)
    - [KMLsToKMLsProperFormat](#kmlstoklmsproperformat)
  - [3. HTML Map Generator Classes](#3-html-map-generator-classes)
    - [InteractiveCategoricalHtmlMap](#interactivecategoricalhtmlmap)
    - [InteractiveContinuousHtmlMap](#interactivecontinuoushtmlmap)
- [Usage Examples](#usage-examples)
  - [StandardGeodataframe](#standardgeodataframe)
  - [KML Converters](#kml-converters)
  - [HTML Map Generators](#html-map-generators)
- [Extensibility](#extensibility)
  - [GeoDataFrameAdapter](#geodataframeadapter)
  - [Adding New Implementations](#adding-new-implementations)
- [Contributing](#contributing)
  - [Development Setup](#development-setup)
  - [Testing](#testing)
  - [Guidelines for Contributors](#guidelines-for-contributors)
  - [Adding New Features](#adding-new-features)
  - [Implementing Support for GeoPolars](#implementing-support-for-geopolars)
- [Documentation](#documentation)
- [License](#license)
- [Author](#author)

## Features

- **Standardize** GeoDataFrames with consistent CRS, valid geometries, and handling of Z-coordinates and polygon holes
- **Convert** between KML files and GeoDataFrames with proper formatting and metadata
- **Visualize** spatial data with interactive HTML maps supporting both categorical and continuous data

## Installation

```bash
pip install geokitten
```

## Dependencies

- `geopandas`: For spatial data operations
- `shapely`: For geometry manipulation
- `numpy`: For numerical operations
- `scipy`: For spatial algorithms
- `folium`: For creating interactive maps
- `branca`: For colormap generation
- `matplotlib`: For color schemes
- Additional dependencies: `os`, `typing`, `pandas`, `fiona`, `xml.etree.ElementTree`, `xml.dom.minidom`

## Quick Start

```python
import geopandas as gpd
from geokitten import (
    StandardGeodataframe,
    KMLsToGeodataframe,
    GeodataframeToKMLs,
    InteractiveCategoricalHtmlMap
)

# Standardize a GeoDataFrame
gdf = gpd.read_file("path/to/shapefile.shp")
standard_gdf = StandardGeodataframe.from_geodataframe(gdf)

# Add surface area in square kilometers
standard_gdf = standard_gdf.get_km2_surface_area()

# Convert KML files to a consolidated GeoDataFrame
kml_converter = KMLsToGeodataframe("path/to/kml_directory")
kml_gdf = kml_converter.consolidate(id_column_name="RegionID")

# Create an interactive map
map_creator = InteractiveCategoricalHtmlMap(
    gdf=standard_gdf,
    color_column="region",
    tooltip_columns=["region", "SURF_A_KM2"]
)
map_creator.create(
    color_scheme="Set1",
    fill_opacity=0.7,
    title="Regional Map",
    output_file="regional_map.html"
)
```

## Module Overview

### 1. StandardGeodataframe Class

The `StandardGeodataframe` class extends the GeoPandas GeoDataFrame to provide standardization and additional functionality for working with geospatial data.

#### Constructor

```python
StandardGeodataframe(gdf_input, crs="EPSG:4326", remove_geni=True, **kwargs)
```

- **gdf_input**: Input GeoDataFrame or file path to be standardized
- **crs**: Coordinate reference system to standardize to (defaults to "EPSG:4326")
- **remove_geni**: Whether to remove "geni" (holes) from polygons (defaults to True)
- **kwargs**: Additional arguments to pass to GeoDataFrame constructor

#### Class Methods

##### `from_file`

```python
StandardGeodataframe.from_file(file_path, crs="EPSG:4326", remove_geni=True, **kwargs)
```

Creates a StandardGeodataframe from a file path.

- **file_path**: Path to the geospatial file (shp, geojson, etc.)
- **crs**: Target CRS (defaults to "EPSG:4326")
- **remove_geni**: Whether to remove polygon holes (defaults to True)
- **kwargs**: Additional arguments for gdp.read_file()

##### `from_geodataframe`

```python
StandardGeodataframe.from_geodataframe(gdf, crs="EPSG:4326", remove_geni=True)
```

Creates a StandardGeodataframe from an existing GeoDataFrame.

- **gdf**: Input GeoDataFrame
- **crs**: Target CRS (defaults to "EPSG:4326")
- **remove_geni**: Whether to remove polygon holes (defaults to True)

#### Methods

##### `get_interior_points`

```python
get_interior_points(geometry_column='geometry', inplace=False, column_name='interior_point')
```

Calculates a point that is guaranteed to be inside each geometry.

- **geometry_column**: Name of the geometry column (defaults to 'geometry')
- **inplace**: Whether to add the interior points to the GeoDataFrame (defaults to False)
- **column_name**: Name of the column to store interior points if inplace=True (defaults to 'interior_point')
- **Returns**: GeoDataFrame with interior points or a GeoSeries of interior points if inplace=False

##### `substract_overlapping_geometries`

```python
substract_overlapping_geometries(column_name, args, remove_geni=True, inplace=False)
```

Subtracts overlapping geometries based on specified criteria.

- **column_name**: Name of the column containing identifiers for geometries
- **args**: Either a tuple ([target_ids], [subtractor_ids]) or a dictionary {target_id: [subtractor_ids]}
- **remove_geni**: Whether to remove holes created by subtraction (defaults to True)
- **inplace**: Whether to modify the GeoDataFrame in place (defaults to False)
- **Returns**: Modified GeoDataFrame if inplace=False

##### `get_m2_surface_area`

```python
get_m2_surface_area(inplace=False, column_name='SURF_A_M2')
```

Calculates the surface area of each geometry in square meters.

- **inplace**: Whether to add the area to the GeoDataFrame (defaults to False)
- **column_name**: Name of the column to store areas if inplace=True (defaults to 'SURF_A_M2')
- **Returns**: GeoDataFrame with area column or a Series of areas if inplace=False

##### `get_km2_surface_area`

```python
get_km2_surface_area(inplace=False, column_name='SURF_A_KM2')
```

Calculates the surface area of each geometry in square kilometers.

- **inplace**: Whether to add the area to the GeoDataFrame (defaults to False)
- **column_name**: Name of the column to store areas if inplace=True (defaults to 'SURF_A_KM2')
- **Returns**: GeoDataFrame with area column or a Series of areas if inplace=False

### 2. KML Converter Classes

#### KMLsToGeodataframe

Converts KML files in a directory into a consolidated GeoDataFrame.

##### Constructor

```python
KMLsToGeodataframe(kml_input_dir)
```

- **kml_input_dir**: Directory containing KML files to process

##### Methods

###### `consolidate`

```python
consolidate(id_column_name='Name', verbose=True)
```

Consolidates all KML files in the directory into a single GeoDataFrame.

- **id_column_name**: Column name for the KML file identifiers (defaults to 'Name')
- **verbose**: Whether to print progress information (defaults to True)
- **Returns**: Consolidated GeoDataFrame containing data from all KML files

#### GeodataframeToKMLs

Converts a GeoDataFrame into individual KML files.

##### Constructor

```python
GeodataframeToKMLs(gdf_input, id_column_name="Name")
```

- **gdf_input**: Input GeoDataFrame or file path
- **id_column_name**: Column name containing identifiers for KML files

##### Methods

###### `kml_metadata_format`

```python
kml_metadata_format(output_dir, verbose=True)
```

Saves each row of the GeoDataFrame as an individual KML file with proper formatting.

- **output_dir**: Directory to save the KML files
- **verbose**: Whether to print progress information (defaults to True)

#### KMLsToKMLsProperFormat

Transforms KML files into a standardized format.

##### Constructor

```python
KMLsToKMLsProperFormat(kml_input_dir, output_dir, id_column_name="Name")
```

- **kml_input_dir**: Directory containing KML files to process
- **output_dir**: Directory to save the transformed KML files
- **id_column_name**: Column name for the KML file identifiers (defaults to 'Name')

##### Methods

###### `transform_format`

```python
transform_format(verbose=True)
```

Transforms KML files into a standardized format and saves them in the output directory.

- **verbose**: Whether to print progress information (defaults to True)

### 3. HTML Map Generator Classes

#### InteractiveCategoricalHtmlMap

Creates categorical interactive choropleth maps.

##### Constructor

```python
InteractiveCategoricalHtmlMap(gdf, color_column, tooltip_columns=None)
```

- **gdf**: GeoDataFrame containing spatial data
- **color_column**: Column name in the GeoDataFrame for categorical coloring
- **tooltip_columns**: List of column names to include in tooltips (defaults to None)

##### Methods

###### `create`

```python
create(color_scheme='tab20', set_custom_colors=None, fill_opacity=0.5, line_opacity=1.0, line_weight=1, location=None, zoom_start=None, title=None, legend_name=None, output_file='categorical_map.html')
```

Generates a categorical choropleth map and saves it as an HTML file.

- **color_scheme**: Color scheme name (e.g., 'tab20', 'Set1') (defaults to 'tab20')
- **set_custom_colors**: List of custom hex color codes (defaults to None)
- **fill_opacity**: Opacity of polygon fill (0-1) (defaults to 0.5)
- **line_opacity**: Opacity of polygon borders (0-1) (defaults to 1.0)
- **line_weight**: Width of polygon borders (defaults to 1)
- **location**: Center coordinates [lat, lon] for map (defaults to None, auto-calculated)
- **zoom_start**: Initial zoom level (defaults to None, auto-calculated)
- **title**: Title for the map (defaults to None)
- **legend_name**: Name for the legend (defaults to None)
- **output_file**: Filename to save the HTML map (defaults to 'categorical_map.html')
- **Returns**: The created folium.Map object

#### InteractiveContinuousHtmlMap

Creates continuous interactive choropleth maps.

##### Constructor

```python
InteractiveContinuousHtmlMap(gdf, color_column, tooltip_columns=None)
```

- **gdf**: GeoDataFrame containing spatial data
- **color_column**: Column name in the GeoDataFrame for continuous coloring
- **tooltip_columns**: List of column names to include in tooltips (defaults to None)

##### Methods

###### `create`

```python
create(color_scheme='viridis', fill_opacity=0.5, line_opacity=1.0, line_weight=1, location=None, zoom_start=None, title=None, legend_name=None, output_file='continuous_map.html')
```

Generates a continuous choropleth map and saves it as an HTML file.

- **color_scheme**: Matplotlib colormap name (e.g., 'viridis', 'coolwarm') (defaults to 'viridis')
- **fill_opacity**: Opacity of polygon fill (0-1) (defaults to 0.5)
- **line_opacity**: Opacity of polygon borders (0-1) (defaults to 1.0)
- **line_weight**: Width of polygon borders (defaults to 1)
- **location**: Center coordinates [lat, lon] for map (defaults to None, auto-calculated)
- **zoom_start**: Initial zoom level (defaults to None, auto-calculated)
- **title**: Title for the map (defaults to None)
- **legend_name**: Name for the legend (defaults to None)
- **output_file**: Filename to save the HTML map (defaults to 'continuous_map.html')
- **Returns**: The created folium.Map object

## Usage Examples

### StandardGeodataframe

```python
from geokitten import StandardGeodataframe
import geopandas as gpd

# Create from file with default settings
gdf = StandardGeodataframe.from_file("path/to/shapefile.shp")

# Create with custom CRS and keeping polygon holes
gdf = StandardGeodataframe.from_file(
    "path/to/shapefile.shp",
    crs="EPSG:3857",
    remove_geni=False
)

# Create from existing GeoDataFrame
gpd_gdf = gpd.read_file("path/to/shapefile.shp")
gdf = StandardGeodataframe.from_geodataframe(gpd_gdf)

# Add interior points as a column
gdf.get_interior_points(inplace=True, column_name="centroid")

# Calculate surface area in square kilometers
gdf_with_area = gdf.get_km2_surface_area(inplace=True)

# Subtract overlapping geometries using tuple arguments
gdf.substract_overlapping_geometries(
    column_name="region_id",
    args=(["target_region"], ["overlapping_region1", "overlapping_region2"]),
    inplace=True
)

# Subtract overlapping geometries using dictionary arguments
subtraction_dict = {
    "target_region1": ["overlapping_region1", "overlapping_region2"],
    "target_region2": ["overlapping_region3"]
}
result = gdf.substract_overlapping_geometries(
    column_name="region_id",
    args=subtraction_dict,
    remove_geni=True,
    inplace=False
)
```

### KML Converters

```python
from geokitten import KMLsToGeodataframe, GeodataframeToKMLs, KMLsToKMLsProperFormat

# Convert KML files to a GeoDataFrame
kml_converter = KMLsToGeodataframe("path/to/kml_directory")
kml_gdf = kml_converter.consolidate(id_column_name="RegionID", verbose=True)

# Convert a GeoDataFrame to KML files
converter = GeodataframeToKMLs(kml_gdf, id_column_name="RegionID")
converter.kml_metadata_format("path/to/output_directory", verbose=True)

# Transform KML files to standardized format
transformer = KMLsToKMLsProperFormat(
    kml_input_dir="path/to/raw_kml_files",
    output_dir="path/to/formatted_kml_files",
    id_column_name="RegionID"
)
transformer.transform_format(verbose=True)
```

### HTML Map Generators

```python
from geokitten import InteractiveCategoricalHtmlMap, InteractiveContinuousHtmlMap

# Create a categorical map
cat_map_creator = InteractiveCategoricalHtmlMap(
    gdf=gdf,
    color_column="region",
    tooltip_columns=["region", "population", "SURF_A_KM2"]
)
cat_map = cat_map_creator.create(
    color_scheme="Set1",
    fill_opacity=0.7,
    line_opacity=0.9,
    line_weight=2,
    title="Regional Map",
    legend_name="Regions",
    output_file="regional_map.html"
)

# Create a continuous map
cont_map_creator = InteractiveContinuousHtmlMap(
    gdf=gdf,
    color_column="population_density",
    tooltip_columns=["region", "population_density", "SURF_A_KM2"]
)
cont_map = cont_map_creator.create(
    color_scheme="coolwarm",
    fill_opacity=0.7,
    title="Population Density Map",
    legend_name="Population Density",
    output_file="density_map.html"
)
```

## Extensibility

GeoKitten is designed with extensibility in mind. While it currently supports GeoPandas for handling geospatial data, the architecture (particularly the `standardize` module) is built to accommodate future extensions like GeoPolar support.

### GeoDataFrameAdapter

The `GeoDataFrameAdapter` class provides a unified interface that can be extended to support other geospatial data frame implementations beyond GeoPandas. This adapter pattern allows the package to evolve to support new spatial dataframe libraries without breaking existing functionality.

```python
from geokitten import StandardGeodataframe

# The current implementation uses GeoPandas under the hood
standard_gdf = StandardGeodataframe.from_file("path/to/shapefile.shp")

# In the future, GeoPolar or other implementations could be supported
# through the same interface without changing user code
```

### Adding New Implementations

To add support for GeoPolar or other spatial dataframe libraries:

1. Extend the `GeoDataFrameAdapter` class with a new implementation
2. Implement the required methods to match the interface
3. Update the `_InputTransformer` class to detect and handle the new implementation

This architecture ensures that as the geospatial Python ecosystem evolves, GeoKitten can adapt without requiring significant changes to user code or the core API.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Development Setup

1. Clone the testing branch repository
2. Install development dependencies:
   ```
    "geopandas>=0.9.0",
    "numpy>=1.23.4",
    "geopandas>=0.14.2",
    "pandas>=2.2.2",
    "shapely>=2.0.2",
    "pytest>=8.3.5",
    "fiona>=1.8.22",
    "folium>=0.14.0",
    "branca>=0.6.0",
    "matplotlib>=3.7.1",
    "pytest>=8.3"
   ```

### Testing

GeoKitten includes a comprehensive test suite located in the `tests/` directory. The tests include:

- `test_standardize.py`: Tests for the StandardGeodataframe class and related functionality
- `test_convert.py`: Tests for KML conversion utilities
- `test_visualize.py`: Tests for HTML map generation

The `tests/tests_files/` directory contains sample data files used for testing.

Before submitting a pull request, please ensure that all tests pass:

```bash
pytest tests/*.py
```

If you're adding new functionality, please add appropriate tests in the `tests/` directory.

### Guidelines for Contributors

1. **Update Tests**: Any new functionality must include corresponding tests
2. **Run Tests**: Ensure all tests pass before submitting a pull request
3. **Documentation**: Follow the existing code style and update documentation
4. **Code Quality**: Follow PEP 8 style guidelines for Python code
5. **Compatibility**: Maintain backward compatibility when possible

### Adding New Features

When adding new functionality, consider the following:

1. Is it consistent with existing APIs?
2. Does it follow the extensibility patterns already established?
3. Is it well-documented with docstrings and example usage?
4. Are there comprehensive tests covering the new functionality?

### Implementing Support for GeoPolars

If you're interested in implementing GeoPolars support:

1. Create a new adapter class that implements the `GeoDataFrameAdapter` interface
2. Update the `_InputTransformer` to detect and handle GeoPolars dataframes
3. Add tests to ensure compatibility with existing functionality
4. Update documentation to include GeoPolars examples

## Documentation

For full documentation, visit [github(GeoKitten)](https://github.com/sape94/GeoKitten).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

- **Sergio A. Pelayo Escalera** - [sergioapelayoe@gmail.com](mailto:sergioapelayoe@gmail.com)
