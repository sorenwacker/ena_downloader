"""
Tests for the ENADownloader class.
"""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from ena_downloader.downloader import ENADownloader


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield tmp_dir


@pytest.fixture
def downloader(temp_dir):
    """Create a downloader instance for tests."""
    return ENADownloader(output_dir=temp_dir, resume=True)


def test_initialization(temp_dir):
    """Test initialization of the downloader."""
    downloader = ENADownloader(output_dir=temp_dir)

    print("Downloader initialized with output directory:", downloader.output_dir)

    # Check if directory structure was created
    assert os.path.exists(os.path.join(temp_dir, "raw", "sequencing"))
    assert os.path.exists(os.path.join(temp_dir, "raw", "metadata"))
    assert os.path.exists(os.path.join(temp_dir, "processed", "aligned"))
    assert os.path.exists(os.path.join(temp_dir, "reference", "genome"))

    # Check if log file was created
    assert os.path.exists(os.path.join(temp_dir, "download_log.txt"))


@patch("requests.get")
def test_search_dataset(mock_get, downloader):
    """Test searching for dataset."""
    # Mock the response
    mock_response = MagicMock()
    mock_response.json.return_value = [
        {
            "study_accession": "PRJEB12345",
            "sample_accession": "SAMEA12345",
            "run_accession": "ERR12345",
            "fastq_ftp": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR12345/ERR12345_1.fastq.gz;ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR12345/ERR12345_2.fastq.gz",
            "fastq_md5": "abcdef1234567890;0987654321fedcba",
        }
    ]
    mock_response.raise_for_status.return_value = None
    mock_get.return_value = mock_response

    # Test the function
    result = downloader.search_dataset("PRJEB12345")

    # Check the result
    assert len(result) == 1
    assert result[0]["study_accession"] == "PRJEB12345"
    assert "fastq_ftp" in result[0]
    assert "fastq_md5" in result[0]

    # Check the API call
    mock_get.assert_called_once()
    args, kwargs = mock_get.call_args
    assert kwargs["params"]["query"] == "study_accession=PRJEB12345"
    assert kwargs["params"]["result"] == "read_run"


@patch("requests.get")
def test_search_dataset_error(mock_get, downloader):
    """Test error handling when searching for dataset."""
    # Mock the response with an error
    mock_response = MagicMock()
    mock_response.raise_for_status.side_effect = Exception("API Error")
    mock_get.return_value = mock_response

    # Test the function
    result = downloader.search_dataset("INVALID")

    # Check the result is an empty list
    assert result == []


def test_verify_md5(temp_dir):
    """Test MD5 verification."""
    # Create a test file
    test_file = os.path.join(temp_dir, "test_file.txt")
    with open(test_file, "w") as f:
        f.write("test content")

    # Calculate the MD5 of "test content"
    # MD5 of "test content" is 9473fdd0d880a43c21b7778d34872157

    downloader = ENADownloader(output_dir=temp_dir)

    # Test with correct MD5
    assert downloader.verify_md5(test_file, "9473fdd0d880a43c21b7778d34872157")

    # Test with incorrect MD5
    assert not downloader.verify_md5(test_file, "incorrect_md5")


@patch("ena_downloader.downloader.ENADownloader.download_file")
@patch("ena_downloader.downloader.ENADownloader.search_dataset")
def test_download_dataset(mock_search, mock_download, downloader):
    """Test downloading a dataset."""
    # Mock the search results
    mock_search.return_value = [
        {
            "study_accession": "PRJEB12345",
            "sample_accession": "SAMEA12345",
            "run_accession": "ERR12345",
            "fastq_ftp": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR12345/ERR12345_1.fastq.gz;ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR12345/ERR12345_2.fastq.gz",
            "fastq_md5": "abcdef1234567890;0987654321fedcba",
            "submitted_ftp": "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR123/ERR12345/ERR12345.txt",
        }
    ]

    # Mock successful downloads
    mock_download.side_effect = ["/tmp/ERR12345_1.fastq.gz", "/tmp/ERR12345_2.fastq.gz", "/tmp/SAMEA12345.txt"]

    # Test the function
    result = downloader.download_dataset("PRJEB12345")

    # Check the result
    assert len(result) == 3
    assert mock_download.call_count == 3


@patch("ena_downloader.downloader.ENADownloader.download_file")
@patch("ena_downloader.downloader.ENADownloader.search_dataset")
def test_download_dataset_with_max_files(mock_search, mock_download, downloader):
    """Test downloading a dataset with max files limit."""
    # Mock the search results
    mock_search.return_value = [
        {
            "study_accession": "PRJEB12345",
            "sample_accession": "SAMEA12345",
            "run_accession": "ERR12345",
            "fastq_ftp": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR12345/ERR12345_1.fastq.gz;ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR123/ERR12345/ERR12345_2.fastq.gz",
            "fastq_md5": "abcdef1234567890;0987654321fedcba",
            "submitted_ftp": "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR123/ERR12345/ERR12345.txt",
        }
    ]

    # Mock successful downloads
    mock_download.return_value = "/tmp/file.fastq.gz"

    # Test the function with max_files=1
    result = downloader.download_dataset("PRJEB12345", max_files=1)

    # Check the result - should only be one file
    assert len(result) == 1
    assert mock_download.call_count == 1
