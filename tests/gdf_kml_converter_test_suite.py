"""
This module contains unit tests and integration tests for the `gdf_kml_converter` module.

It tests the functionality of the `KMLsToGeodataframe`, `_KMLsProperFormat`, `GeodataframeToKMLs`, 
and `KMLsToKMLsProperFormat` classes, ensuring proper handling of KML files, GeoDataFrame 
transformations, and KML file generation.

Author: Sergio A. Pelayo Escalera
Version: 0.1.0
Created: 2025-04-09
Last Updated: 2025-04-09
"""

from geokitten.gdf_kml_converter import (
    KMLsToGeodataframe,
    _KMLsProperFormat,
    GeodataframeToKMLs,
    KMLsToKMLsProperFormat
)
__version__ = "0.1.0"
__author__ = "Sergio A. Pelayo Escalera (sergio.pelayo@nielseniq.com)"
__collabs__ = [""]
__created__ = "2025-04-09"
__last_updated__ = "2025-04-09"

import os
import pytest
import geopandas as gpd
from unittest.mock import patch, MagicMock, mock_open
import xml.etree.ElementTree as ET
import sys
from pathlib import Path

# Add the src directory to the Python path
sys.path.append(str(Path(__file__).parent.parent))


class TestKMLsToGeodataframe:
    """
    Unit tests for the `KMLsToGeodataframe` class.

    These tests validate the functionality of reading KML files from a directory,
    consolidating them into a GeoDataFrame, and handling invalid inputs.
    """

    @pytest.fixture
    def kml_dir(self):
        """Fixture that returns the path to test KML files"""
        return "./tests/tests_files/inputs/gdf_kml_converter/KMLsToGeodataframe"

    @pytest.fixture
    def validation_gdf(self):
        """Fixture that returns the validation GeoDataFrame"""
        validation_file = "./tests/tests_files/outputs/gdf_kml_converter_KMLsToGeodataframe_consolidate_test_validation_file.shp"
        return gpd.read_file(validation_file)

    def test_initialization(self, kml_dir):
        """
        Test successful initialization of `KMLsToGeodataframe`.

        Ensures that the class is initialized correctly with a valid directory.
        """
        converter = KMLsToGeodataframe(kml_dir)
        assert converter.is_valid is None  # Should not raise an error
        assert len(converter.kml_files) > 0
        assert converter.kml_dir == kml_dir

    def test_kml_dir_validator_valid(self, kml_dir):
        """
        Test `_kml_dir_validator` with a valid directory.

        Ensures that no exception is raised for a valid directory.
        """
        converter = KMLsToGeodataframe(kml_dir)
        assert converter._kml_dir_validator(kml_dir) is None

    def test_kml_dir_validator_invalid(self):
        """
        Test `_kml_dir_validator` with an invalid directory.

        Ensures that a `ValueError` is raised for a nonexistent directory.
        """
        with pytest.raises(ValueError, match="Directory .* does not exist"):
            KMLsToGeodataframe("./nonexistent_directory")

    def test_kml_files_validator_valid(self, kml_dir):
        """
        Test `_kml_files_validator` with a directory containing KML files.

        Ensures that the correct list of KML files is returned.
        """
        converter = KMLsToGeodataframe(kml_dir)
        kml_files = converter._kml_files_validator(kml_dir)
        assert len(kml_files) > 0
        assert all(f.endswith('.kml') for f in kml_files)

    def test_kml_files_validator_no_files(self, tmp_path):
        """
        Test `_kml_files_validator` with a directory containing no KML files.

        Ensures that a `ValueError` is raised if no KML files are found.
        """
        with pytest.raises(ValueError, match="No KML files found in directory"):
            KMLsToGeodataframe(str(tmp_path))

    def test_safely_read_kml_remove_geni_false(self, kml_dir):
        """
        Test `_safely_read_kml` with a valid KML file.
        Setting remove_geni to False (default).

        Ensures that the KML file is correctly read into a GeoDataFrame.
        """
        converter = KMLsToGeodataframe(kml_dir)
        kml_file = converter.kml_files[0]
        kml_gdf = converter._safely_read_kml(
            kml_file=kml_file,
            kml_file_name="test_name",
            verbose=False)
        assert isinstance(kml_gdf, gpd.GeoDataFrame)
        assert "Name" in kml_gdf.columns
        assert "geometry" in kml_gdf.columns
        assert len(kml_gdf.columns) == 2

    def test_safely_read_kml_remove_geni_true(self, kml_dir):
        """
        Test `_safely_read_kml` with a valid KML file.
        Setting remove_geni to True.

        Ensures that the KML file is correctly read into a GeoDataFrame.
        """
        converter = KMLsToGeodataframe(kml_dir)
        kml_file = converter.kml_files[0]
        kml_gdf = converter._safely_read_kml(
            kml_file, "test_name", remove_geni=True, verbose=False)
        assert isinstance(kml_gdf, gpd.GeoDataFrame)
        assert "Name" in kml_gdf.columns
        assert "geometry" in kml_gdf.columns
        assert len(kml_gdf.columns) == 2

    @patch('geokitten.gdf_kml_converter.StandardGeodataframe')
    def test_safely_read_kml_exception(self, mock_std_gdf, kml_dir):
        """
        Test `_safely_read_kml` with an invalid KML file.

        Ensures that `None` is returned and an error message is printed.
        """
        mock_std_gdf.side_effect = Exception("Test exception")
        converter = KMLsToGeodataframe(kml_dir)
        result = converter._safely_read_kml(
            "test.kml", "test_name", verbose=False)
        assert result is None

    def test_gdf_list_append(self, kml_dir):
        """
        Test `_gdf_list_append` method.

        Ensures that a GeoDataFrame is correctly appended to the list of GeoDataFrames.
        """
        converter = KMLsToGeodataframe(kml_dir)
        kml_file = converter.kml_files[0]
        kml_gdf = converter._safely_read_kml(
            kml_file, "test_name", verbose=False)
        # Test with valid GDF
        kml_gdf_list = []
        result_list = converter._gdf_list_append(
            kml_gdf_list, kml_gdf, verbose=False)
        assert len(result_list) == 1
        assert result_list[0] is kml_gdf
        # Test with None GDF
        kml_gdf_list = []
        result_list = converter._gdf_list_append(
            kml_gdf_list, None, verbose=False)
        assert len(result_list) == 0

    def test_concat_kml_gdf_list(self, kml_dir):
        """
        Test `_concat_kml_gdf_list` method.

        Ensures that a list of GeoDataFrames is correctly consolidated into a single GeoDataFrame.
        """
        converter = KMLsToGeodataframe(kml_dir)
        kml_gdf_list = []
        # Create a list of GDFs
        # Use just the first two files for testing
        for kml_file in converter.kml_files[:2]:
            kml_gdf = converter._safely_read_kml(
                kml_file, kml_file.replace('.kml', ''), verbose=False)
            kml_gdf_list = converter._gdf_list_append(
                kml_gdf_list, kml_gdf, verbose=False)

        result_gdf = converter._concat_kml_gdf_list(
            kml_gdf_list, verbose=False)
        assert isinstance(result_gdf, gpd.GeoDataFrame)
        assert len(result_gdf) == sum(len(gdf) for gdf in kml_gdf_list)
        assert result_gdf.crs == kml_gdf_list[0].crs

    def test_consolidate_remove_geni_false(self, kml_dir, validation_gdf):
        """
        Test `consolidate` method against a validation file.
        Setting remove_geni to False (default).

        Ensures that all KML files are properly read and consolidated into a single GeoDataFrame.
        """
        converter = KMLsToGeodataframe(kml_dir)
        result_gdf = converter.consolidate('Id', verbose=False)
        # Check basic properties
        assert isinstance(result_gdf, gpd.GeoDataFrame)
        assert "Id" in result_gdf.columns
        assert "geometry" in result_gdf.columns
        # Compare with validation GDF
        # We should check that all the KML files were properly read and consolidated
        # This comparison depends on the actual data, so we'll check some basic properties
        assert result_gdf.crs == validation_gdf.crs
        assert len(result_gdf) == len(validation_gdf)
        # Check that all names from validation file exist in result
        validation_names = set(validation_gdf["Id"].tolist())
        result_names = set(result_gdf["Id"].tolist())
        assert validation_names == result_names

    def test_consolidate_remove_geni_true(self, kml_dir, validation_gdf):
        """
        Test `consolidate` method against a validation file.
        Setting remove_geni to False (default).

        Ensures that all KML files are properly read and consolidated into a single GeoDataFrame.
        """
        converter = KMLsToGeodataframe(kml_dir)
        result_gdf = converter.consolidate(
            'Id', remove_geni=True, verbose=False)
        assert isinstance(result_gdf, gpd.GeoDataFrame)
        assert "Id" in result_gdf.columns
        assert "geometry" in result_gdf.columns
        assert result_gdf.crs == validation_gdf.crs
        assert len(result_gdf) == len(validation_gdf)
        validation_names = set(validation_gdf["Id"].tolist())
        result_names = set(result_gdf["Id"].tolist())
        assert validation_names == result_names

    def test_consolidate_custom_id_column(self, kml_dir):
        """
        Test `consolidate` method with a custom ID column name.

        Ensures that the custom column name is correctly used in the consolidated GeoDataFrame.
        """
        converter = KMLsToGeodataframe(kml_dir)
        result_gdf = converter.consolidate(
            id_column_name="CustomID", verbose=False)
        assert "CustomID" in result_gdf.columns


class Test_KMLsProperFormat:
    """
    Unit tests for the `_KMLsProperFormat` class.

    These tests validate the functionality of formatting GeoDataFrames as KML files,
    generating KML structures, and saving them to the output directory.
    """

    @pytest.fixture
    def sample_gdf(self):
        """
        Fixture that returns a sample GeoDataFrame.

        Provides a GeoDataFrame for testing `_KMLsProperFormat`.
        """
        validation_file = "./tests/tests_files/outputs/gdf_kml_converter_KMLsToGeodataframe_consolidate_test_validation_file.shp"
        return gpd.read_file(validation_file)

    @pytest.fixture
    def output_dir(self, tmp_path):
        """
        Fixture that returns a temporary directory for output.

        Provides a temporary directory for saving KML files during tests.
        """
        output_path = tmp_path / "output"
        output_path.mkdir()
        return str(output_path)

    def test_initialization(self, sample_gdf, output_dir):
        """
        Test initialization of `_KMLsProperFormat`.

        Ensures that the class is initialized correctly with a GeoDataFrame and output directory.
        """
        formatter = _KMLsProperFormat(sample_gdf, output_dir, verbose=False)
        assert formatter.kml_gdf is sample_gdf
        assert formatter.output_dir == output_dir
        assert 'Id' in formatter.kml_gdf.columns

    def test_validate_or_create_output_dir_existing(self, sample_gdf, output_dir):
        """
        Test `_validate_or_create_output_dir` with an existing directory.

        Ensures that the existing directory is correctly validated.
        """
        formatter = _KMLsProperFormat(sample_gdf, output_dir, verbose=False)
        result_dir = formatter._validate_or_create_output_dir(
            output_dir, verbose=False)
        assert result_dir == output_dir
        assert os.path.exists(output_dir)

    def test_validate_or_create_output_dir_new(self, sample_gdf, tmp_path):
        """
        Test `_validate_or_create_output_dir` with a new directory.

        Ensures that a new directory is created if it does not exist.
        """
        new_dir = str(tmp_path / "new_output")
        formatter = _KMLsProperFormat(sample_gdf, new_dir, verbose=False)
        assert os.path.exists(new_dir)

    @patch.object(_KMLsProperFormat, '_generate_kml_structure')
    @patch.object(_KMLsProperFormat, '_add_polygons')
    @patch.object(_KMLsProperFormat, '_format_kml')
    @patch('builtins.open', new_callable=mock_open)
    def test_save_kml(self, mock_file, mock_format_kml, mock_add_polygons,
                      mock_generate_kml, sample_gdf, output_dir):
        """
        Test `save_kml` method.

        Ensures that the KML file is correctly saved to the output directory.
        """
        # Setup mocks
        mock_kml = MagicMock()
        mock_generate_kml.return_value = mock_kml
        mock_document = MagicMock()
        mock_kml.find.return_value = mock_document
        mock_format_kml.return_value = "formatted_kml_content"
        # Test just one row of the sample_gdf to avoid complex mocking
        test_gdf = sample_gdf.iloc[:1].copy()
        formatter = _KMLsProperFormat(test_gdf, output_dir, verbose=False)
        formatter.save_kml(verbose=False)
        # Verify method calls
        mock_generate_kml.assert_called_once()
        mock_kml.find.assert_called_once_with('Document')
        mock_add_polygons.assert_called_once_with(mock_document)
        mock_format_kml.assert_called_once_with(mock_kml)
        # Check file operations
        mock_file.assert_called_once()
        mock_file().write.assert_called_once_with("formatted_kml_content")

    def test_generate_kml_structure(self, sample_gdf, output_dir):
        """
        Test `_generate_kml_structure` method.

        Ensures that the base KML structure is correctly generated.
        """
        test_gdf = sample_gdf.iloc[:1].copy()  # Use just one row
        formatter = _KMLsProperFormat(
            test_gdf, output_dir, id_column_name='Id', verbose=False)
        # Need to patch 'gdf' attribute manually since it's used in _generate_kml_structure
        formatter.gdf = test_gdf
        kml = formatter._generate_kml_structure()
        # Basic structure checks
        assert kml.tag == 'kml'
        assert kml.attrib['xmlns'] == 'http://earth.google.com/kml/2.2'
        document = kml.find('Document')
        assert document is not None
        # Check if the document name is set correctly
        document_name = document.find('name')
        assert document_name is not None
        assert document_name.text == test_gdf['Name'].iloc[0]
        # Check if Document is set to open
        open_element = document.find('open')
        assert open_element is not None
        assert open_element.text == '1'

    def test_add_styles(self, sample_gdf, output_dir):
        """
        Test `_add_styles` method.

        Ensures that styles are correctly added to the KML document.
        """
        formatter = _KMLsProperFormat(sample_gdf, output_dir, verbose=False)
        # Create a document element to add styles to
        document = ET.Element('Document')
        formatter._add_styles(document)
        # Check if block style was added
        block_style = document.find("Style[@id='for_block_styling']")
        assert block_style is not None
        # Check if sub block style was added
        sub_block_style = document.find("Style[@id='for_sub_block_styling']")
        assert sub_block_style is not None
        # Check if style maps were added
        style_map_block = document.find(
            "StyleMap[@id='sty_for_block_styling']")
        assert style_map_block is not None
        style_map_sub_block = document.find(
            "StyleMap[@id='sty_for_sub_block_styling']")
        assert style_map_sub_block is not None

    def test_format_kml(self, sample_gdf, output_dir):
        """
        Test `_format_kml` method.

        Ensures that the KML structure is correctly formatted as a string.
        """
        formatter = _KMLsProperFormat(sample_gdf, output_dir, verbose=False)
        # Create a simple KML structure for testing
        kml = ET.Element('kml')
        document = ET.SubElement(kml, 'Document')
        ET.SubElement(document, 'name').text = 'Test Document'
        formatted_kml = formatter._format_kml(kml)
        # Check if the formatting is correct
        assert '<?xml version="1.0" encoding="UTF-8"?>' in formatted_kml
        assert '<kml xmlns="http://earth.google.com/kml/2.2">' in formatted_kml
        assert '<Document>' in formatted_kml
        assert '<name>Test Document</name>' in formatted_kml
        assert '</Document>' in formatted_kml
        assert '</kml>' in formatted_kml


class TestGeodataframeToKMLs:
    """
    Unit tests for the `GeodataframeToKMLs` class.

    These tests validate the functionality of converting a GeoDataFrame into individual
    KML files and handling invalid inputs.
    """

    @pytest.fixture
    def input_gdf_path(self):
        """Fixture that returns the path to test GeoDataFrame"""
        return "./tests/tests_files/inputs/gdf_kml_converter_GeodataframeToKMLS_test_file.shp"

    @pytest.fixture
    def validation_dir(self):
        """Fixture that returns the path to validation KML files"""
        return "./tests/tests_files/outputs/gdf_kml_converter/GeodataframeToKMLs"

    @pytest.fixture
    def sample_gdf(self, input_gdf_path):
        """Fixture that returns the sample GeoDataFrame"""
        return gpd.read_file(input_gdf_path)

    def test_initialization_with_file(self, input_gdf_path):
        """
        Test initialization with a file path.

        Ensures that the class is correctly initialized with a valid file path.
        """
        with patch('geokitten.gdf_kml_converter.StandardGeodataframe') as mock_std_gdf:
            mock_gdf = MagicMock()
            mock_gdf.columns = ['Id']
            mock_gdf['Id'].is_unique = True
            mock_std_gdf.return_value = mock_gdf
            converter = GeodataframeToKMLs(
                input_gdf_path, id_column_name='Id', remove_geni=True)
            assert converter.id_column_name == 'Id'
            assert converter.gdf is mock_gdf

    def test_initialization_with_gdf(self, sample_gdf):
        """
        Test initialization with a GeoDataFrame.

        Ensures that the class is correctly initialized with a valid GeoDataFrame.
        """
        with patch('geokitten.gdf_kml_converter.StandardGeodataframe') as mock_std_gdf:
            mock_std_gdf.return_value = sample_gdf
            converter = GeodataframeToKMLs(sample_gdf, id_column_name='Id')
            assert converter.id_column_name == 'Id'
            assert converter.gdf is sample_gdf

    def test_validate_id_column_name_existence_valid(self, sample_gdf):
        """
        Test `_validate_id_column_name_existence` with a valid column.

        Ensures that no exception is raised for an existing column.
        """
        with patch('geokitten.gdf_kml_converter.StandardGeodataframe') as mock_std_gdf:
            mock_std_gdf.return_value = sample_gdf
            # Add a Name column to the sample_gdf if it doesn't exist
            if 'Id' not in sample_gdf.columns:
                sample_gdf['Id'] = [f"Id_{i}" for i in range(len(sample_gdf))]
            converter = GeodataframeToKMLs(sample_gdf, id_column_name='Id')
            result = converter._validate_id_column_name_existence('Id')
            assert result is None

    def test_validate_id_column_name_existence_invalid(self, sample_gdf):
        """
        Test `_validate_id_column_name_existence` with an invalid column.

        Ensures that a `ValueError` is raised for a nonexistent column.
        """
        with patch('geokitten.gdf_kml_converter.StandardGeodataframe') as mock_std_gdf:
            mock_std_gdf.return_value = sample_gdf
            converter = GeodataframeToKMLs(sample_gdf, id_column_name='Id')
            with pytest.raises(ValueError, match="Column 'NonExistentColumn' not found in GeoDataFrame"):
                converter._validate_id_column_name_existence(
                    'NonExistentColumn')

    def test_validate_uniqueness_valid(self, sample_gdf):
        """
        Test `_validate_uniqueness` with unique values.

        Ensures that no exception is raised for a column with unique values.
        """
        with patch('geokitten.gdf_kml_converter.StandardGeodataframe') as mock_std_gdf:
            # Ensure Name column has unique values
            if 'Name' not in sample_gdf.columns:
                sample_gdf['Name'] = [
                    f"Name_{i}" for i in range(len(sample_gdf))]
            else:
                sample_gdf['Name'] = [
                    f"{sample_gdf['Name'].iloc[i]}_{i}" for i in range(len(sample_gdf))]
            mock_std_gdf.return_value = sample_gdf
            # This should not raise an error
            converter = GeodataframeToKMLs(sample_gdf)
            assert converter.id_column_name == 'Name'

    def test_validate_uniqueness_invalid(self, sample_gdf):
        """
        Test `_validate_uniqueness` with non-unique values.

        Ensures that a `ValueError` is raised for a column with duplicate values.
        """
        with patch('geokitten.gdf_kml_converter.StandardGeodataframe') as mock_std_gdf:
            # Make Name column have duplicate values
            if 'Name' not in sample_gdf.columns:
                sample_gdf['Name'] = ['Same_Name'] * len(sample_gdf)
            else:
                sample_gdf['Name'] = ['Same_Name'] * len(sample_gdf)
            mock_std_gdf.return_value = sample_gdf
            with pytest.raises(ValueError, match="Column Name must have unique values"):
                GeodataframeToKMLs(sample_gdf)

    @patch.object(GeodataframeToKMLs, '_update_id_column')
    @patch('geokitten.gdf_kml_converter._KMLsProperFormat')
    def test_kml_metadata_format(self, mock_kml_proper_format, mock_catch_warnings, sample_gdf):
        """
        Test `_kml_metadata_format` method.

        Ensures that each row of the GeoDataFrame is saved as an individual KML file.
        """
        with patch('geokitten.gdf_kml_converter.StandardGeodataframe') as mock_std_gdf:
            # Ensure Name column has unique values
            if 'Name' not in sample_gdf.columns:
                sample_gdf['Name'] = [
                    f"Name_{i}" for i in range(len(sample_gdf))]
            else:
                sample_gdf['Name'] = [
                    f"{sample_gdf['Name'].iloc[i]}_{i}" for i in range(len(sample_gdf))]
            mock_std_gdf.return_value = sample_gdf
            converter = GeodataframeToKMLs(sample_gdf)
            # Mock the necessary components
            mock_kml_instance = MagicMock()
            mock_kml_proper_format.return_value = mock_kml_instance
            # Call the method
            converter.kml_metadata_format("test_output_dir", verbose=False)
            # Verify calls
            assert mock_catch_warnings.call_count == len(sample_gdf)
            assert mock_kml_proper_format.call_count == len(sample_gdf)
            assert mock_kml_instance.save_kml.call_count == len(sample_gdf)


class TestKMLsToKMLsProperFormat:
    """
    Unit tests for the `KMLsToKMLsProperFormat` class.

    These tests validate the functionality of transforming KML files into a standardized
    format and saving them as properly formatted KML files.
    """

    @pytest.fixture
    def input_dir(self):
        """Fixture that returns the path to test KML files"""
        return "./tests/tests_files/inputs/gdf_kml_converter/KMLsToKMLsProperFormat"

    @pytest.fixture
    def validation_dir(self):
        """Fixture that returns the path to validation KML files"""
        return "./tests/tests_files/outputs/gdf_kml_converter/KMLsToKMLsProperFormat"

    @pytest.fixture
    def output_dir(self, tmp_path):
        """Fixture that returns a temporary directory for output"""
        output_path = tmp_path / "output"
        output_path.mkdir()
        return str(output_path)

    def test_initialization(self, input_dir, output_dir):
        """
        Test initialization of `KMLsToKMLsProperFormat`.

        Ensures that the class is correctly initialized with valid input and output directories.
        """
        with patch('geokitten.gdf_kml_converter.KMLsToGeodataframe') as mock_kmls_to_gdf:
            mock_kmls_to_gdf.return_value._kml_files_validator.return_value = [
                'test1.kml', 'test2.kml']
            converter = KMLsToKMLsProperFormat(input_dir, output_dir)
            assert converter.input_dir == input_dir
            assert converter.output_dir == output_dir
            assert converter.id_col == 'Name'
            assert converter.kml_files == ['test1.kml', 'test2.kml']

    def test_initialization_default_output(self, input_dir):
        """
        Test initialization with the default output directory.

        Ensures that the output directory defaults to the input directory if not specified.
        """
        with patch('geokitten.gdf_kml_converter.KMLsToGeodataframe') as mock_kmls_to_gdf:
            mock_kmls_to_gdf.return_value._kml_files_validator.return_value = [
                'test1.kml', 'test2.kml']
            converter = KMLsToKMLsProperFormat(input_dir)
            assert converter.input_dir == input_dir
            assert converter.output_dir == input_dir

    def test_path_validation_valid(self, input_dir):
        """
        Test `_path_validation` with a valid path.

        Ensures that the directory path is correctly validated.
        """
        converter = KMLsToKMLsProperFormat(input_dir, input_dir)
        # Should not raise an error
        converter._path_validation(input_dir)

    def test_path_validation_invalid(self):
        """
        Test `_path_validation` with an invalid path.

        Ensures that a `ValueError` is raised for a nonexistent directory.
        """
        with pytest.raises(ValueError, match="Directory .* does not exist"):
            KMLsToKMLsProperFormat("./nonexistent_directory")

    def test_validate_or_create_output_dir_existing(self, input_dir):
        """
        Test `_validate_or_create_output_dir` with an existing directory.

        Ensures that the existing directory is correctly validated.
        """
        converter = KMLsToKMLsProperFormat(input_dir, input_dir)
        result_dir = converter._validate_or_create_output_dir(
            input_dir, verbose=False)
        assert result_dir == input_dir

    def test_validate_or_create_output_dir_new(self, input_dir, tmp_path):
        """
        Test `_validate_or_create_output_dir` with a new directory.

        Ensures that a new directory is created if it does not exist.
        """
        new_dir = str(tmp_path / "new_output")
        with patch('geokitten.gdf_kml_converter.KMLsToGeodataframe'):
            converter = KMLsToKMLsProperFormat(input_dir, new_dir)
            assert os.path.exists(new_dir)

    @patch('geokitten.gdf_kml_converter._KMLsProperFormat')
    def test_transform_format(self, mock_kml_proper_format, input_dir, output_dir):
        """
        Test `transform_format` method.

        Ensures that KML files are correctly transformed into a standardized format
        and saved to the output directory.
        """
        # Setup mocks
        mock_kmls_to_gdf = MagicMock()
        mock_kml_gdf = MagicMock()
        mock_kmls_to_gdf._safely_read_kml.return_value = mock_kml_gdf

        # Create test instance with mocked KMLsToGeodataframe
        with patch('geokitten.gdf_kml_converter.KMLsToGeodataframe', return_value=mock_kmls_to_gdf):
            mock_kmls_to_gdf._kml_files_validator.return_value = [
                'test1.kml', 'test2.kml']
            converter = KMLsToKMLsProperFormat(
                input_dir, output_dir, id_column_name='Id')
            converter.transform_format(verbose=False)
            # Verify calls
            assert mock_kmls_to_gdf._safely_read_kml.call_count == 2
            # Check that each call uses 'Id' as the id_column_name
            for call_args in mock_kmls_to_gdf._safely_read_kml.call_args_list:
                assert call_args[0][2] == 'Id'
            assert mock_kml_proper_format.call_count == 2
            # Check that each call uses 'Id' as the id_column_name
            for call_args in mock_kml_proper_format.call_args_list:
                assert call_args[0][2] == 'Id'
            # Verify that save_kml was called on each _KMLsProperFormat instance
            mock_proper_format_instance = mock_kml_proper_format.return_value
            assert mock_proper_format_instance.save_kml.call_count == 2
