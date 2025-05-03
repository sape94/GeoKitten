# GeoDataFrame Standardization Module Documentation

## Overview

The `gdf_standardization` module provides classes and methods for standardizing and processing GeoDataFrames. It includes functionality for transforming input GeoDataFrames or file paths into standardized formats, removing Z-coordinates, handling invalid geometries, and processing polygons to remove "geni" (holes). The module supports the GeoPandas framework.

- **Author**: Sergio A. Pelayo Escalera
- **Email**: sergioapelayoe@gmail.com
- **Version:** 0.1.0
- **Created**: 2025-04-09
- **Last Modified**: 2025-04-09

## Table of Contents

- [Dependencies](#dependencies)
- [Classes](#classes)
  - [GeoDataFrameAdapter](#geodataframeadapter)
  - [StandardGeodataframe](#standardgeodataframe)
  - [\_InputTransformer](#_inputtransformer)
  - [\_GeniRemover](#_geniremover)
- [Functions](#functions)
- [Usage Examples](#usage-examples)
- [API Reference](#api-reference)

## Dependencies

The module depends on the following Python packages:

- `geopandas` (for spatial data operations)
- `shapely` (for geometry manipulation)
- `numpy` (for numerical operations)
- `scipy` (for spatial algorithms)
- Standard library packages: `os`, `typing`

```python
import os
import numpy as np
from typing import Union, List, Tuple, Optional, Any
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, LinearRing, GeometryCollection
from scipy.spatial import Voronoi
```

## Classes

### GeoDataFrameAdapter

Adapter class to provide a unified interface for GeoPandas GeoDataFrames.

```python
adapter = GeoDataFrameAdapter(my_geodataframe)
```

#### Key Features

- Provides a consistent interface for GeoPandas GeoDataFrames
- Handles CRS transformations
- Supports common GeoDataFrame operations like copying, indexing, and iteration
- Offers conversion to/from native GeoDataFrame types

### StandardGeodataframe

The main class for standardizing and processing GeoDataFrames. It inherits from `gpd.GeoDataFrame` and provides additional methods for validating geometries, columns, and performing operations such as subtracting overlapping geometries.

```python
# Create from file
gdf = StandardGeodataframe.from_file("path/to/shapefile.shp")

# Create from existing GeoDataFrame
gdf = StandardGeodataframe.from_geodataframe(my_geodataframe)

# Create directly
gdf = StandardGeodataframe(my_geodataframe)
```

#### Key Features

- Standardizes GeoDataFrames to EPSG:4326 CRS (configurable)
- Removes Z-coordinates from geometries
- Processes polygons to remove holes ("geni")
- Validates geometries and fixes invalid ones
- Provides methods to subtract overlapping geometries
- Calculates surface areas in square meters or square kilometers
- Gets interior points for geometries

### \_InputTransformer

Internal class for transforming and standardizing input GeoDataFrames or file paths.

#### Key Features

- Validates input GeoDataFrames or file paths
- Standardizes CRS
- Removes Z-coordinates
- Optionally removes "geni" (holes) from polygons

### \_GeniRemover

Internal class for removing "geni" (holes) from a Polygon geometry.

#### Key Features

- Processes polygons to merge holes with the exterior ring
- Identifies and handles the nearest hole to the exterior ring
- Creates "cuts" from the exterior ring to holes

## Functions

### extend_geoseries_with_interior_point

Extends the GeoPandas GeoSeries class with an 'interior_point' property.

```python
# This is called automatically when the module is imported
extend_geoseries_with_interior_point()

# After this, any GeoSeries can use the property
interior_points = my_geoseries.interior_point
```

### \_get_interior_point

Get a point that is guaranteed to be inside a geometry.

```python
point = _get_interior_point(polygon)
```

## Usage Examples

### Basic Usage

```python
import geopandas as gpd
from gdf_standardization import StandardGeodataframe

# Create a StandardGeodataframe from a file
gdf = StandardGeodataframe.from_file("path/to/shapefile.shp")

# Create a StandardGeodataframe from an existing GeoDataFrame
gpd_gdf = gpd.read_file("path/to/shapefile.shp")
gdf = StandardGeodataframe.from_geodataframe(gpd_gdf)

# Calculate surface area in square kilometers
gdf_with_area = gdf.get_km2_surface_area()

# Get interior points for each geometry
interior_points = gdf.get_interior_points()

# Add interior points as a column
gdf.get_interior_points(inplace=True, column_name="centroid")
```

### Subtracting Overlapping Geometries

```python
# Example with tuple arguments
gdf.substract_overlapping_geometries(
    column_name="region_id",
    args=(["target_region"], ["overlapping_region1", "overlapping_region2"]),
    inplace=True
)

# Example with dictionary arguments
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

### Using the GeoDataFrameAdapter

```python
from gdf_standardization import GeoDataFrameAdapter

# Create an adapter for a GeoDataFrame
adapter = GeoDataFrameAdapter(my_geodataframe)

# Transform to a different CRS
transformed = adapter.to_crs("EPSG:3857")

# Access the native GeoDataFrame
native_gdf = adapter.to_native()
```

## API Reference

### GeoDataFrameAdapter

```python
GeoDataFrameAdapter(geodataframe)
```

- `to_crs(crs, inplace=False)`: Transform the GeoDataFrame to a new CRS
- `crs`: Property to get or set the CRS
- `copy()`: Create a deep copy of the GeoDataFrame
- `columns`: Property to get the column names
- `__getitem__(key)`: Get items from the GeoDataFrame
- `iterrows()`: Iterate over rows as (index, Series) pairs
- `apply(func, axis=0)`: Apply a function along an axis
- `to_native()`: Return the native GeoDataFrame object
- `from_file(path)`: Class method to create an adapter from a file

### StandardGeodataframe

```python
StandardGeodataframe(gdf_input, crs="EPSG:4326", remove_geni=True, **kwargs)
```

- `from_file(file_path, crs="EPSG:4326", remove_geni=True, **kwargs)`: Create from a file
- `from_geodataframe(gdf, crs="EPSG:4326", remove_geni=True)`: Create from an existing GeoDataFrame
- `get_interior_points(geometry_column='geometry', inplace=False, column_name='interior_point')`: Get interior points
- `substract_overlapping_geometries(column_name, args, remove_geni=True, inplace=False)`: Subtract overlapping geometries
- `get_m2_surface_area(inplace=False, column_name='SURF_A_M2')`: Calculate surface area in square meters
- `get_km2_surface_area(inplace=False, column_name='SURF_A_KM2')`: Calculate surface area in square kilometers

### \_InputTransformer (Internal)

```python
_InputTransformer(gdf_input, crs="EPSG:4326", remove_geni=True)
```

- `set_standard_crs()`: Standardize the CRS of the GeoDataFrame
- `apply_z_coord_and_geni_removal(gdf)`: Apply Z-coordinate and hole removal

### \_GeniRemover (Internal)

```python
_GeniRemover(geom)
```

- `trnsf_pol_all_geni()`: Transform the Polygon by removing all holes

### Extension Functions

- `extend_geoseries_with_interior_point()`: Extends GeoSeries with 'interior_point' property
- `_get_interior_point(geometry)`: Get a point inside a geometry
- `_get_interior_point_property(self)`: Property implementation for interior points
