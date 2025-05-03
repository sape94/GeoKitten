# KML to GeoDataFrame Converter

A Python module for converting KML files to GeoDataFrames and vice versa, as well as transforming KML files into a standardized format.

- **Author**: Sergio A. Pelayo Escalera
- **Email**: sergioapelayoe@gmail.com
- **Version**: 0.1.0
- **Created**: 2025-04-09
- **Last Modified**: 2025-04-09

## Overview

This module provides classes and methods for handling KML (Keyhole Markup Language) files commonly used in geographic information systems. It facilitates:

- Reading KML files and converting them to GeoDataFrames
- Consolidating multiple KML files into a single GeoDataFrame
- Converting GeoDataFrames back to KML files with proper formatting
- Transforming existing KML files into a standardized format

## Requirements

The module has the following dependencies:

- pandas
- geopandas
- fiona
- xml.etree.ElementTree
- xml.dom.minidom

## Classes

### KMLsToGeodataframe

A class to convert KML files in a directory into a consolidated GeoDataFrame.

#### Initialization

```python
from src.gdf_kml_converter import KMLsToGeodataframe

# Initialize with directory containing KML files
kml_converter = KMLsToGeodataframe("path/to/kml_directory")
```

#### Methods

##### `consolidate(id_column_name='Name', verbose=True)`

Consolidates all KML files in the directory into a single GeoDataFrame.

**Parameters**:

- `id_column_name` (str, optional): The column name for the KML file name. Defaults to 'Name'.
- `verbose` (bool, optional): Whether to print verbose output. Defaults to True.

**Returns**:

- `gpd.GeoDataFrame`: The consolidated GeoDataFrame.

**Example**:

```python
# With default column name
consolidated_gdf = kml_converter.consolidate()
print(consolidated_gdf.head())

# With custom column name
consolidated_gdf = kml_converter.consolidate(id_column_name="CustomID")
print(consolidated_gdf.head())
```

### GeodataframeToKMLs

A class to convert a GeoDataFrame into individual KML files.

#### Initialization

```python
from src.gdf_kml_converter import GeodataframeToKMLs

# Initialize with a GeoDataFrame
converter = GeodataframeToKMLs(gdf, id_column_name="Name")

# Or initialize with a file path
converter = GeodataframeToKMLs("path/to/shapefile.shp", id_column_name="Name")
```

#### Methods

##### `kml_metadata_format(output_dir, verbose=True)`

Saves each row of the GeoDataFrame as an individual KML file with proper formatting.

**Parameters**:

- `output_dir` (str): The directory to save the KML files.
- `verbose` (bool, optional): Whether to print verbose output. Defaults to True.

**Example**:

```python
import geopandas as gpd
from shapely.geometry import Point

# Create a sample GeoDataFrame
geometry = [Point(0, 0), Point(1, 1)]
gdf = gpd.GeoDataFrame({'Name': ['Point1', 'Point2'], 'geometry': geometry}, crs="EPSG:4326")

# Convert to KML files
converter = GeodataframeToKMLs(gdf, id_column_name="Name")
converter.kml_metadata_format(output_dir="path/to/output_directory", verbose=True)
```

### KMLsToKMLsProperFormat

A class to transform KML files into a standardized format.

#### Initialization

```python
from src.gdf_kml_converter import KMLsToKMLsProperFormat

# Initialize with input and output directories
transformer = KMLsToKMLsProperFormat(
    kml_input_dir="path/to/input_directory",
    output_dir="path/to/output_directory",
    id_column_name="CustomID"
)
```

#### Methods

##### `transform_format(verbose=True)`

Transforms KML files into a standardized format and saves them in the output directory.

**Parameters**:

- `verbose` (bool, optional): Whether to print verbose output. Defaults to True.

**Example**:

```python
transformer = KMLsToKMLsProperFormat(
    kml_input_dir="path/to/input_directory",
    output_dir="path/to/output_directory",
    id_column_name="CustomID"
)
transformer.transform_format(verbose=True)
```

## Helper Classes

### \_KMLsProperFormat

An internal class used by `GeodataframeToKMLs` and `KMLsToKMLsProperFormat` to format and save GeoDataFrames as KML files with proper styling.

This class handles:

- Generating the KML structure with proper XML formatting
- Adding styles for polygons (including line and fill styles)
- Creating placemark elements for each geometry

## Usage Examples

### Converting KML Files to a GeoDataFrame

```python
from src.gdf_kml_converter import KMLsToGeodataframe

# Initialize the converter
kml_converter = KMLsToGeodataframe("path/to/kml_files")

# Consolidate KML files into a single GeoDataFrame
consolidated_gdf = kml_converter.consolidate(id_column_name="RegionID")

# Use the GeoDataFrame for analysis or visualization
print(consolidated_gdf.head())
```

### Converting a GeoDataFrame to KML Files

```python
from src.gdf_kml_converter import GeodataframeToKMLs
import geopandas as gpd

# Load a GeoDataFrame from a shapefile
gdf = gpd.read_file("path/to/shapefile.shp")

# Initialize the converter
converter = GeodataframeToKMLs(gdf, id_column_name="RegionName")

# Save each row as a KML file
converter.kml_metadata_format("path/to/output_directory")
```

### Standardizing KML Files

```python
from src.gdf_kml_converter import KMLsToKMLsProperFormat

# Initialize the transformer
transformer = KMLsToKMLsProperFormat(
    kml_input_dir="path/to/raw_kml_files",
    output_dir="path/to/formatted_kml_files"
)

# Transform and save the KML files
transformer.transform_format()
```

## Notes

- All KML geometries are expected to be in WGS 84 (EPSG:4326) coordinate system.
- The module handles both Polygon and MultiPolygon geometries.
- KML files are formatted with standard Google Earth styling, including proper line colors and transparency.
