"""
Geokitten: A Python package for geospatial data processing and visualization.
"""
__version__ = '0.1.0'
__author__ = "Sergio A. Pelayo Escalera (sergioapelayoe@gmail.com)"
__collabs__ = [""]
__created__ = "2025-04-09"
__last_updated__ = "2025-04-09"

# Import main functions to expose at package level
from .gdf_standardization import StandardGeodataframe
from .gdf_kml_converter import (
    KMLsToGeodataframe,
    GeodataframeToKMLs,
    KMLsToKMLsProperFormat
)
from .html_generator import (
    InteractiveCategoricalHtmlMap,
    InteractiveContinuousHtmlMap
)

# Define what should be imported with "from geokit import *"
__all__ = [
    'StandardGeodataframe',
    'KMLsToGeodataframe',
    'GeodataframeToKMLs',
    'KMLsToKMLsProperFormat',
    'InteractiveCategoricalHtmlMap',
    'InteractiveContinuousHtmlMap'
]
