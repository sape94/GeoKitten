"""
This module provides classes and methods for converting KML files to GeoDataFrames
and vice versa, as well as for transforming KML files into a standardized format.

It includes functionality for reading KML files, consolidating them into GeoDataFrames,
and saving GeoDataFrames back into KML files with proper formatting.

Author: Sergio A. Pelayo Escalera
Email: sergioapelayoe@gmail.com
Created: 2025-04-09
Last-modified: 2025-04-09
Version: 0.1.0
"""

from geokitten.gdf_standardization import StandardGeodataframe
__version__ = "0.1.0"
__author__ = "Sergio A. Pelayo Escalera (sergioapelayoe@gmail.com)"
__collabs__ = [""]
__created__ = "2025-04-09"
__last_updated__ = "2025-04-09"

import os
import pandas as pd
import geopandas as gpd
from typing import Union
from xml.dom import minidom
import xml.etree.ElementTree as ET
import fiona

fiona.supported_drivers['KML'] = 'rw'


class KMLsToGeodataframe:
    """
    A class to convert KML files in a directory into a consolidated GeoDataFrame.

    This class validates the input directory, reads KML files, and consolidates
    them into a single GeoDataFrame.

    Attributes:
        is_valid (bool): Indicates whether the input directory is valid.
        kml_files (list): A list of KML file names in the directory.
        kml_dir (str): The input directory containing KML files.
    """

    def __init__(self,
                 kml_dir: str):
        """
        Initialize the KMLsToGeodataframe instance.

        This constructor sets up the instance by validating the input directory 
        and retrieving the list of KML files within it.

        Args:
            kml_dir (str): The directory containing KML files.

        Raises:
            ValueError: If the directory does not exist or contains no KML files.

        Examples:
            Example: Initializing the KMLsToGeodataframe instance
            -----------------------------------------------------
            >>> from src.gdf_kml_converter import KMLsToGeodataframe
            >>> kml_converter = KMLsToGeodataframe("path/to/kml_directory")
            >>> print(kml_converter.kml_files)
        """
        self.is_valid = self._kml_dir_validator(kml_dir)
        self.kml_files = self._kml_files_validator(kml_dir)
        self.kml_dir = kml_dir

    def _kml_dir_validator(self,
                           kml_dir: str) -> None:
        """
        Validate the existence of the input directory.

        Args:
            kml_dir (str): The directory to validate.

        Raises:
            ValueError: If the directory does not exist.
        """
        if not os.path.exists(kml_dir):
            raise ValueError(f"Directory {kml_dir} does not exist.")

    def _kml_files_validator(self,
                             kml_dir: str) -> list:
        """
        Validate and retrieve KML files in the directory.

        Args:
            kml_dir (str): The directory to search for KML files.

        Returns:
            list: A list of KML file names in the directory.

        Raises:
            ValueError: If no KML files are found in the directory.
        """
        kml_files = [f for f in os.listdir(kml_dir) if f.endswith('.kml')]
        if len(kml_files) == 0:
            raise ValueError(f"No KML files found in directory {kml_dir}.")
        return kml_files

    def _safely_read_kml(self,
                         kml_file: str,
                         kml_file_name: str,
                         id_column_name: str = 'Name',
                         remove_geni: bool = False,
                         verbose: bool = True) -> Union[gpd.GeoDataFrame, None]:
        """
        Safely read a KML file into a GeoDataFrame.

        Args:
            kml_file (str): The KML file to read.
            kml_file_name (str): The name of the KML file (without extension).
            id_column_name (str, optional): The column name for the KML file name. Defaults to 'Name'.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Returns:
            Union[gpd.GeoDataFrame, None]: The GeoDataFrame or None if an error occurs.
        """
        try:
            kml_gdf = gpd.read_file(
                f'{self.kml_dir}/{kml_file}', driver='KML')
            kml_gdf = StandardGeodataframe(kml_gdf,
                                           remove_geni=remove_geni)
            kml_gdf[f'{id_column_name}'] = kml_file_name
            kml_gdf = kml_gdf[[f'{id_column_name}', 'geometry']]
            return kml_gdf
        except Exception as e:
            if verbose:
                print(f"Error reading KML file {kml_file}: {e}")
            return None

    def _gdf_list_append(self,
                         kml_gdf_list: list,
                         kml_gdf: gpd.GeoDataFrame,
                         verbose: bool = True) -> list:
        """
        Append a GeoDataFrame to the list of GeoDataFrames.

        Args:
            kml_gdf_list (list): The list of GeoDataFrames.
            kml_gdf (gpd.GeoDataFrame): The GeoDataFrame to append.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Returns:
            list: The updated list of GeoDataFrames.
        """
        if kml_gdf is not None:
            kml_gdf_list.append(kml_gdf)
            if verbose:
                print(f"\t✓ Successfully read.")
        return kml_gdf_list

    def _concat_kml_gdf_list(self,
                             kml_gdf_list: list,
                             verbose: bool = True) -> gpd.GeoDataFrame:
        """
        Consolidate a list of GeoDataFrames into a single GeoDataFrame.

        Args:
            kml_gdf_list (list): The list of GeoDataFrames.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Returns:
            gpd.GeoDataFrame: The consolidated GeoDataFrame.
        """
        consolidated_gdf = pd.concat(kml_gdf_list, ignore_index=True)
        consolidated_gdf = gpd.GeoDataFrame(consolidated_gdf,
                                            crs=kml_gdf_list[0].crs)
        if verbose:
            print(f"✓✓✓ KML files consolidated successfully.")
        return consolidated_gdf

    def consolidate(self,
                    id_column_name: str = 'Name',
                    remove_geni: bool = False,
                    verbose: bool = True) -> gpd.GeoDataFrame:
        """
        Consolidate all KML files in the directory into a single GeoDataFrame.

        This method reads all KML files in the specified directory, processes them into
        GeoDataFrames, and consolidates them into a single GeoDataFrame. Each KML file's
        name is stored in the specified column.

        Args:
            id_column_name (str, optional): The column name for the KML file name. Defaults to 'Name'.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Returns:
            gpd.GeoDataFrame: The consolidated GeoDataFrame.

        Examples:
            Example 1: Consolidating KML files with default column name
            -----------------------------------------------------------
            >>> from src.gdf_kml_converter import KMLsToGeodataframe
            >>> kml_converter = KMLsToGeodataframe("path/to/kml_directory")
            >>> consolidated_gdf = kml_converter.consolidate()
            >>> print(consolidated_gdf.head())

            Example 2: Consolidating KML files with a custom column name
            ------------------------------------------------------------
            >>> from src.gdf_kml_converter import KMLsToGeodataframe
            >>> kml_converter = KMLsToGeodataframe("path/to/kml_directory")
            >>> consolidated_gdf = kml_converter.consolidate(id_column_name="CustomID")
            >>> print(consolidated_gdf.head())
        """
        kml_gdf_list = []
        for kml_file in self.kml_files:
            kml_file_name = kml_file.replace('.kml', '')
            if verbose:
                print(f'Processing: {kml_file_name}')
            kml_gdf = self._safely_read_kml(kml_file,
                                            kml_file_name,
                                            id_column_name,
                                            remove_geni,
                                            verbose)
            kml_gdf_list = self._gdf_list_append(
                kml_gdf_list, kml_gdf, verbose)
        consolidated_gdf = self._concat_kml_gdf_list(kml_gdf_list, verbose)
        return consolidated_gdf


class _KMLsProperFormat:
    """
    A class to format and save GeoDataFrames as KML files.

    This class provides methods to generate KML structures, add styles, and save
    GeoDataFrames as properly formatted KML files.

    Attributes:
        kml_gdf (gpd.GeoDataFrame): The GeoDataFrame to save as a KML file.
        output_dir (str): The directory to save the KML file.
    """

    def __init__(self,
                 kml_gdf: gpd.GeoDataFrame,
                 output_dir: str,
                 id_column_name: str = 'Name',
                 verbose: bool = True) -> None:
        """
        Initialize the _KMLsProperFormat instance.

        Args:
            kml_gdf (gpd.GeoDataFrame): The GeoDataFrame to save as a KML file.
            output_dir (str): The directory to save the KML file.
            id_column_name (str, optional): The column name for the KML file name. Defaults to 'Name'.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.
        """
        self.kml_gdf = kml_gdf
        self.output_dir = self._validate_or_create_output_dir(output_dir,
                                                              verbose)
        self.kml_gdf.rename(
            columns={f'{id_column_name}': 'Name'}, inplace=True)

    def _validate_or_create_output_dir(self,
                                       output_dir: str,
                                       verbose: bool = True) -> str:
        """
        Validate or create the output directory.

        Args:
            output_dir (str): The directory to validate or create.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Returns:
            str: The validated or newly created output directory path.
        """
        if not os.path.exists(output_dir):
            if verbose:
                print(f"Creating output directory: {output_dir}")
            os.makedirs(output_dir, exist_ok=True)
        return output_dir

    def _generate_kml_structure(self) -> ET.Element:
        """
        Generate the base KML structure.

        This method creates the root KML element and a Document element, including
        the document name and styles.

        Returns:
            ET.Element: The root KML element.
        """
        kml = ET.Element('kml')
        kml.set('xmlns', 'http://earth.google.com/kml/2.2')
        document = ET.SubElement(kml, 'Document')
        document_name = self.kml_gdf['Name'].iloc[0] if 'Name' in self.kml_gdf.columns else 'Document'
        ET.SubElement(document, 'name').text = document_name
        ET.SubElement(document, 'open').text = '1'
        self._add_styles(document)
        return kml

    def _add_block_style(self,
                         document: ET.Element) -> None:
        """
        Add block style to the KML document.

        Args:
            document (ET.Element): The Document element to add the style to.
        """
        style_block = ET.SubElement(document, 'Style', id='for_block_styling')
        line_style_block = ET.SubElement(style_block, 'LineStyle')
        ET.SubElement(line_style_block, 'color').text = 'ff0000ff'
        ET.SubElement(line_style_block, 'width').text = '2'
        poly_style_block = ET.SubElement(style_block, 'PolyStyle')
        ET.SubElement(poly_style_block, 'fill').text = '0'

    def _add_sub_block_style(self,
                             document: ET.Element) -> None:
        """
        Add sub-block style to the KML document.

        Args:
            document (ET.Element): The Document element to add the style to.
        """
        style_sub_block = ET.SubElement(
            document, 'Style', id='for_sub_block_styling')
        line_style_sub_block = ET.SubElement(style_sub_block, 'LineStyle')
        ET.SubElement(line_style_sub_block, 'color').text = 'ff0000ff'
        ET.SubElement(line_style_sub_block, 'width').text = '2'
        poly_style_sub_block = ET.SubElement(style_sub_block, 'PolyStyle')
        ET.SubElement(poly_style_sub_block, 'fill').text = '0'

    def _add_style_maps(self,
                        document: ET.Element) -> None:
        """
        Add style maps to the KML document.

        Args:
            document (ET.Element): The Document element to add the style maps to.
        """
        style_map_block = ET.SubElement(
            document, 'StyleMap', id='sty_for_block_styling')
        self._add_style_map_pairs(style_map_block, '#for_block_styling')
        style_map_sub_block = ET.SubElement(
            document, 'StyleMap', id='sty_for_sub_block_styling')
        self._add_style_map_pairs(
            style_map_sub_block, '#for_sub_block_styling')

    def _add_styles(self,
                    document: ET.Element) -> None:
        """
        Add all styles to the KML document.

        Args:
            document (ET.Element): The Document element to add the styles to.
        """
        self._add_block_style(document)
        self._add_sub_block_style(document)
        self._add_style_maps(document)

    def _add_style_map_pairs(self,
                             style_map: ET.Element, style_url: str) -> None:
        """
        Add style map pairs to a style map element.

        Args:
            style_map (ET.Element): The StyleMap element to add pairs to.
            style_url (str): The URL of the style to reference.
        """
        pair_normal = ET.SubElement(style_map, 'Pair')
        ET.SubElement(pair_normal, 'key').text = 'normal'
        ET.SubElement(pair_normal, 'styleUrl').text = style_url
        pair_highlight = ET.SubElement(style_map, 'Pair')
        ET.SubElement(pair_highlight, 'key').text = 'highlight'
        ET.SubElement(pair_highlight, 'styleUrl').text = style_url

    def _add_polygons(self,
                      document: ET.Element) -> None:
        """
        Add polygons from the GeoDataFrame to the KML document.

        Args:
            document (ET.Element): The Document element to add polygons to.
        """
        folder = ET.SubElement(document, 'Folder')
        document_name = self.kml_gdf['Name'].iloc[0] if 'Name' in self.kml_gdf.columns else 'Document'
        ET.SubElement(folder, 'name').text = document_name
        for idx, row in self.kml_gdf.iterrows():
            self._add_placemark(folder, row)

    def _add_placemark(self,
                       folder: ET.Element,
                       row: gpd.GeoSeries) -> None:
        """
        Add a placemark for a row in the GeoDataFrame.

        Args:
            folder (ET.Element): The Folder element to add the placemark to.
            row (gpd.GeoSeries): The row from the GeoDataFrame.
        """
        placemark = ET.SubElement(folder, 'Placemark')
        placemark_name = row['Name'] if 'Name' in row else 'Placemark'
        ET.SubElement(placemark, 'name').text = placemark_name
        ET.SubElement(placemark, 'styleUrl').text = '#sty_for_block_styling'
        if row['geometry'].geom_type == 'Polygon':
            polygons = [row['geometry']]
        elif row['geometry'].geom_type == 'MultiPolygon':
            # Access the .geoms attribute
            polygons = list(row['geometry'].geoms)
        else:
            return  # Skip unsupported geometry types
        for polygon_geom in polygons:
            self._add_polygon(placemark, polygon_geom)

    def _add_polygon(self,
                     placemark: ET.Element,
                     polygon_geom) -> None:
        """
        Add a polygon geometry to a placemark.

        Args:
            placemark (ET.Element): The Placemark element to add the polygon to.
            polygon_geom: The polygon geometry to add.
        """
        polygon = ET.SubElement(placemark, 'Polygon')
        outer = ET.SubElement(polygon, 'outerBoundaryIs')
        linear_ring = ET.SubElement(outer, 'LinearRing')
        ET.SubElement(linear_ring, 'tessellate').text = '1'
        coords = '\n'.join([f'{x},{y},0.0000' for x, y,
                           *_ in polygon_geom.exterior.coords])
        ET.SubElement(linear_ring, 'coordinates').text = f'\n{coords}\n'

    def _format_kml(self,
                    kml: ET.Element) -> str:
        """
        Format the KML structure into a string.

        Args:
            kml (ET.Element): The root KML element.

        Returns:
            str: The formatted KML string.
        """
        xmlstr = minidom.parseString(ET.tostring(kml)).toprettyxml(indent="  ")
        start = xmlstr.find('<Document>')
        end = xmlstr.find('</Document>') + len('</Document>')
        document_content = xmlstr[start:end]
        return f'''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.2">
{document_content}
</kml>'''

    def save_kml(self,
                 verbose: bool = True) -> None:
        """
        Save the GeoDataFrame as a KML file.

        Args:
            verbose (bool, optional): Whether to print verbose output. Defaults to True.
        """
        kml = self._generate_kml_structure()
        document = kml.find('Document')
        self._add_polygons(document)
        formatted_kml = self._format_kml(kml)
        document_name = self.kml_gdf['Name'].iloc[0] if 'Name' in self.kml_gdf.columns else 'Document'
        document_name = document_name.replace('.kml', '')
        output_path = os.path.join(self.output_dir, f'{document_name}.kml')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(formatted_kml)
        if verbose:
            print(f"KML file saved to: {output_path}")


class GeodataframeToKMLs:
    """
    A class to convert a GeoDataFrame into individual KML files.

    This class validates the GeoDataFrame and ensures that each row is saved as
    a separate KML file with proper formatting.

    Attributes:
        gdf (gpd.GeoDataFrame): The input GeoDataFrame.
        id_column_name (str): The column name used to identify individual KML files.
    """

    def __init__(self,
                 gdf_input: Union[gpd.GeoDataFrame, str],
                 id_column_name: str = 'Name',
                 remove_geni: bool = False) -> None:
        """
        Initialize the GeodataframeToKMLs instance.

        This constructor sets up the instance by validating the input GeoDataFrame or file path 
        and ensuring that the specified column name exists and contains unique values.

        Args:
            gdf_input (Union[gpd.GeoDataFrame, str]): The input GeoDataFrame or file path.
            id_column_name (str, optional): The column name used to identify individual KML files. 
                Defaults to 'Name'.
            remove_geni (bool, optional): Whether to remove the 'geni' column from the GeoDataFrame.
                Defaults to False.

        Raises:
            ValueError: If the `id_column_name` does not exist in the GeoDataFrame or contains 
                non-unique values.

        Examples:
            Example 1: Initializing with a GeoDataFrame
            -------------------------------------------
            >>> from src.gdf_kml_converter import GeodataframeToKMLs
            >>> import geopandas as gpd
            >>> from shapely.geometry import Point
            >>> geometry = [Point(0, 0), Point(1, 1)]
            >>> gdf = gpd.GeoDataFrame({'Name': ['Point1', 'Point2'], 'geometry': geometry}, crs="EPSG:4326")
            >>> converter = GeodataframeToKMLs(gdf, id_column_name="Name", remove_geni=False)
            >>> print(converter.gdf)

            Example 2: Initializing with a file path
            ----------------------------------------
            >>> from src.gdf_kml_converter import GeodataframeToKMLs
            >>> converter = GeodataframeToKMLs("path/to/shapefile.shp", id_column_name="Name", remove_geni=False)
            >>> print(converter.gdf)
        """
        self.gdf = StandardGeodataframe(gdf_input, remove_geni=remove_geni)
        self.id_column_name = self._validate_uniqueness(id_column_name)

    def _validate_id_column_name_existence(self,
                                           id_column_name: str) -> str:
        """
        Validate the existence of the ID column in the GeoDataFrame.

        Args:
            id_column_name (str): The column name to validate.

        Raises:
            ValueError: If the column does not exist in the GeoDataFrame.
        """
        if id_column_name not in self.gdf.columns:
            raise ValueError(
                f"Column '{id_column_name}' not found in GeoDataFrame.")

    def _validate_uniqueness(self,
                             id_column_name: str) -> None:
        """
        Validate that the ID column contains unique values.

        Args:
            id_column_name (str): The column name to validate.

        Returns:
            str: The validated column name.

        Raises:
            ValueError: If the column contains non-unique values.
        """
        self._validate_id_column_name_existence(id_column_name)
        if not self.gdf[f'{id_column_name}'].is_unique:
            raise ValueError(
                f"Column {id_column_name} must have unique values.")
        return id_column_name

    def _update_id_column(self,
                          kml_file_name: str) -> None:
        """
        Update the ID column for a specific row.

        Args:
            kml_file_name (str): The name of the KML file to assign to the first row.
        """
        self.gdf.at[self.gdf.index[0], self.id_column_name] = kml_file_name

    def kml_metadata_format(self,
                            output_dir: str,
                            verbose: bool = True) -> str:
        """
        Save each row of the GeoDataFrame as an individual KML file.

        This method iterates through each row of the GeoDataFrame and saves it as an 
        individual KML file in the specified output directory. The KML files are 
        formatted using the `_KMLsProperFormat` class.

        Args:
            output_dir (str): The directory to save the KML files.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Examples:
            Example 1: Saving KML files from a GeoDataFrame
            -----------------------------------------------
            >>> from src.gdf_kml_converter import GeodataframeToKMLs
            >>> import geopandas as gpd
            >>> from shapely.geometry import Point
            >>> geometry = [Point(0, 0), Point(1, 1)]
            >>> gdf = gpd.GeoDataFrame({'Name': ['Point1', 'Point2'], 'geometry': geometry}, crs="EPSG:4326")
            >>> converter = GeodataframeToKMLs(gdf, id_column_name="Name")
            >>> converter.kml_metadata_format(output_dir="path/to/output_directory", verbose=True)
        """
        for idx, row in self.gdf.iterrows():
            self._update_id_column(row[f'{self.id_column_name}'])
            kml = _KMLsProperFormat(self.gdf.loc[[idx]],
                                    output_dir,
                                    self.id_column_name,
                                    verbose)
            kml.save_kml(verbose)


class KMLsToKMLsProperFormat:
    """
    A class to transform KML files into a standardized format.

    This class reads KML files from a directory, processes them into a standardized
    format, and saves them back as properly formatted KML files.

    Attributes:
        input_dir (str): The input directory containing KML files.
        output_dir (str): The output directory to save formatted KML files.
        id_col (str): The column name used to identify individual KML files.
    """

    def __init__(self,
                 kml_input_dir: str,
                 output_dir: str = None,
                 id_column_name: str = 'Name') -> None:
        """
        Initialize the KMLsToKMLsProperFormat instance.

        This constructor sets up the instance by validating the input directory, 
        creating or validating the output directory, and initializing the 
        `KMLsToGeodataframe` instance to handle KML files.

        Args:
            kml_input_dir (str): The input directory containing KML files.
            output_dir (str, optional): The output directory to save formatted KML files. 
                Defaults to the input directory.
            id_column_name (str, optional): The column name used to identify individual 
                KML files. Defaults to 'Name'.

        Raises:
            ValueError: If the input directory does not exist.

        Examples:
            Example: Initializing the KMLsToKMLsProperFormat instance
            ---------------------------------------------------------
            >>> from src.gdf_kml_converter import KMLsToKMLsProperFormat
            >>> transformer = KMLsToKMLsProperFormat(
            ...     kml_input_dir="path/to/input_directory",
            ...     output_dir="path/to/output_directory",
            ...     id_column_name="CustomID"
            ... )
            >>> print(transformer.input_dir)
            >>> print(transformer.output_dir)
        """
        self.input_dir = self._path_validation(kml_input_dir)
        if output_dir is None:
            output_dir = kml_input_dir
        self.output_dir = self._validate_or_create_output_dir(output_dir)
        self.id_col = id_column_name
        self.kmls_gdf_instance = KMLsToGeodataframe(self.input_dir)
        self.kml_files = self.kmls_gdf_instance._kml_files_validator(
            self.input_dir)

    def _path_validation(self,
                         dir_path: str) -> None:
        """
        Validate the existence of a directory.

        Args:
            dir_path (str): The directory path to validate.

        Returns:
            str: The validated directory path.

        Raises:
            ValueError: If the directory does not exist.
        """
        if not os.path.exists(dir_path):
            raise ValueError(f"Directory {dir} does not exist.")
        return dir_path

    def _validate_or_create_output_dir(self,
                                       output_dir: str,
                                       verbose: bool = True) -> str:
        """
        Validate or create the output directory.

        Args:
            output_dir (str): The directory to validate or create.
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Returns:
            str: The validated or newly created output directory path.
        """
        if not os.path.exists(output_dir):
            if verbose:
                print(f"Creating output directory: {output_dir}")
            os.makedirs(output_dir, exist_ok=True)
        return output_dir

    def transform_format(self,
                         remove_geni: bool = False,
                         verbose: bool = True) -> None:
        """
        Transform KML files into a standardized format and save them.

        This method reads KML files from the input directory, processes them into
        a standardized format using the `_KMLsProperFormat` class, and saves them
        back as properly formatted KML files in the output directory.

        Args:
            verbose (bool, optional): Whether to print verbose output. Defaults to True.

        Examples:
            Example: Transforming and saving KML files
            ------------------------------------------
            >>> from src.gdf_kml_converter import KMLsToKMLsProperFormat
            >>> transformer = KMLsToKMLsProperFormat(
            ...     kml_input_dir="path/to/input_directory",
            ...     output_dir="path/to/output_directory",
            ...     id_column_name="CustomID"
            ... )
            >>> transformer.transform_format(verbose=True)
        """
        for kml_file in self.kml_files:
            kml_gdf = self.kmls_gdf_instance._safely_read_kml(kml_file,
                                                              kml_file,
                                                              self.id_col,
                                                              remove_geni,
                                                              verbose)
            self.proper_format = _KMLsProperFormat(kml_gdf,
                                                   self.output_dir,
                                                   self.id_col,
                                                   verbose)
            self.proper_format.save_kml(verbose)
