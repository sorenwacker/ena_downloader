"""
Core implementation of the ENADownloader class.
"""

import ftplib
import hashlib
import logging
import os
import random
from urllib.parse import urlparse

import requests

logger = logging.getLogger(__name__)


class ENADownloader:
    """
    A utility for downloading datasets from the European Nucleotide Archive (ENA).

    Manages downloading of genomic datasets with support for various accession types.
    Creates a folder structure to organize the downloaded files:

    data/
    ├── raw/
    │   ├── sequencing/
    │   │   ├── sample1_R1.fastq.gz
    │   │   ├── sample1_R2.fastq.gz
    │   │   ├── sample2_R1.fastq.gz
    │   │   ├── sample2_R2.fastq.gz
    │   │   └── ...
    │   └── metadata/
    │       ├── sample_metadata.csv
    │       └── ...
    ├── processed/
    │   ├── aligned/
    │   │   ├── sample1.bam
    │   │   ├── sample2.bam
    │   │   └── ...
    │   ├── filtered/
    │   │   ├── sample1_filtered.bam
    │   │   ├── sample2_filtered.bam
    │   │   └── ...
    │   ├── qc/
    │   │   ├── fastqc/
    │   │   │   ├── sample1_R1_fastqc.html
    │   │   │   ├── sample1_R2_fastqc.html
    │   │   │   ├── sample2_R1_fastqc.html
    │   │   │   ├── sample2_R2_fastqc.html
    │   │   │   └── ...
    │   │   └── multiqc/
    │   │       └── multiqc_report.html
    │   ├── counts/
    │   │   ├── sample1_counts.txt
    │   │   ├── sample2_counts.txt
    │   │   └── ...
    │   └── assemblies/
    │       ├── sample1/
    │       │   ├── sample1_contigs.fasta
    │       │   ├── sample1_scaffolds.fasta
    │       │   └── ...
    │       ├── sample2/
    │       │   ├── sample2_contigs.fasta
    │       │   ├── sample2_scaffolds.fasta
    │       │   └── ...
    │       └── ...
    └── reference/
        ├── genome/
        │   ├── genome.fa
        │   └── genome.fa.fai
        └── annotation/
            ├── genes.gtf
            └── ...
    """

    ENA_SEARCH_URL: str = "https://www.ebi.ac.uk/ena/portal/api/search"
    ENA_DOWNLOAD_URL: str = "https://www.ebi.ac.uk/ena/portal/api/download"

    def __init__(self, output_dir: str = "data", resume: bool = True, verbose: bool = False):
        """
        Initialize the ENA downloader.

        Args:
            output_dir: Directory to save downloaded files.
                        Creates the folder structure under this directory.
            resume: Whether to resume downloads of existing files.
                   If True, will skip files that already exist and appear complete.
            verbose: Whether to enable verbose logging.
        """
        self.output_dir = output_dir
        self.resume = resume
        self.verbose = verbose

        # Set up logging if not already configured
        if not logging.getLogger().handlers:
            log_level = logging.DEBUG if verbose else logging.INFO
            logging.basicConfig(level=log_level, format="%(asctime)s - %(levelname)s - %(message)s")

        self.create_folder_structure()

        # Create a file to keep track of downloaded files with their md5 checksums
        self.download_log_file = os.path.join(self.output_dir, "download_log.txt")
        self._initialize_download_log()

    def _initialize_download_log(self):
        """Initialize the download log file if it doesn't exist."""
        if not os.path.exists(self.download_log_file):
            # Create an empty file
            with open(self.download_log_file, "w") as f:
                f.write("# Download log file with MD5 checksums\n")
                f.write("# format: file_path,md5_checksum,status\n")

    def _get_downloaded_files(self) -> dict[str, dict[str, str]]:
        """
        Get a dictionary of already downloaded files and their MD5 checksums.

        Returns:
            Dictionary mapping file paths to their MD5 checksums and status.
        """
        downloaded_files = {}

        if os.path.exists(self.download_log_file):
            with open(self.download_log_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split(",")
                    if len(parts) >= 3:
                        file_path, md5_checksum, status = parts[:3]
                        downloaded_files[file_path] = {"md5": md5_checksum, "status": status}

        return downloaded_files

    def _update_download_log(self, file_path: str, md5_checksum: str, status: str = "complete"):
        """
        Update the download log file with a new entry.

        Args:
            file_path: Path to the downloaded file
            md5_checksum: MD5 checksum of the file
            status: Status of the download ("complete", "incomplete", etc.)
        """
        with open(self.download_log_file, "a") as f:
            f.write(f"{file_path},{md5_checksum},{status}\n")

    def create_folder_structure(self):
        """
        Create the folder structure for organizing the downloaded files.
        """
        # Raw data directories
        os.makedirs(os.path.join(self.output_dir, "raw", "sequencing"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "raw", "metadata"), exist_ok=True)

        # Processed data directories
        os.makedirs(os.path.join(self.output_dir, "processed", "aligned"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "processed", "filtered"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "processed", "qc", "fastqc"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "processed", "qc", "multiqc"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "processed", "counts"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "processed", "assemblies"), exist_ok=True)

        # Reference data directories
        os.makedirs(os.path.join(self.output_dir, "reference", "genome"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "reference", "annotation"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "reference", "metadata"), exist_ok=True)

    def search_dataset(self, accession: str, result_type: str = "read_run") -> list[dict[str, str | int]]:
        """
        Search for dataset details in ENA.

        Args:
            accession: Study or project accession (e.g., PRJNA123456).
            result_type: Type of result to retrieve.

        Returns:
            List of dataset entries.
        """
        params = {
            "result": result_type,
            "query": f"study_accession={accession}",
            "fields": "study_accession,sample_accession,run_accession,fastq_ftp,submitted_ftp,fastq_md5",
            "format": "JSON",
        }

        try:
            logger.info("Searching for dataset: %s", accession)
            response = requests.get(self.ENA_SEARCH_URL, params=params)
            response.raise_for_status()

            dataset_info = response.json()
            logger.info("Found %d entries for %s", len(dataset_info), accession)

            return dataset_info
        except requests.RequestException as e:
            logger.error("Error searching dataset: %s", e)
            logger.error("Request URL: %s", e.request.url if e.request else "N/A")
            logger.error("Response: %s", e.response.text if e.response else "N/A")
            return []

    def download_file_ftp(
        self, url: str, filename: str | None = None, subdir: str = "sequencing", expected_md5: str | None = None
    ) -> str | None:
        """
        Download a file using FTP protocol.

        Args:
            url: FTP URL of the file to download.
            filename: Optional custom filename.
            subdir: Subdirectory under 'raw' to save the file.
            expected_md5: Expected MD5 checksum for verification.

        Returns:
            Path to the downloaded file or None if download fails.
        """
        try:
            # Parse the FTP URL
            parsed_url = urlparse(url)
            host = parsed_url.hostname
            remote_path = parsed_url.path

            # Extract filename if not provided
            if not filename:
                filename = os.path.basename(remote_path)

            full_path = os.path.join(self.output_dir, "raw", subdir, filename)

            # Check if file already exists and if we should resume
            downloaded_files = self._get_downloaded_files()
            if self.resume and full_path in downloaded_files:
                file_info = downloaded_files[full_path]
                if file_info["status"] == "complete":
                    if expected_md5 and file_info["md5"] == expected_md5:
                        logger.info("File %s already exists with matching MD5, skipping download", filename)
                        return full_path
                    if os.path.exists(full_path):
                        # Verify the existing file's MD5
                        if expected_md5 and self.verify_md5(full_path, expected_md5):
                            logger.info("File %s already exists with matching MD5, skipping download", filename)
                            return full_path
                        logger.warning("File %s exists but MD5 does not match, redownloading", filename)
                    else:
                        logger.warning("File %s is in log but not found on disk, redownloading", filename)

            # Check if file exists on disk but not in the log
            if self.resume and os.path.exists(full_path) and expected_md5:
                if self.verify_md5(full_path, expected_md5):
                    logger.info("File %s already exists with matching MD5, skipping download", filename)
                    self._update_download_log(full_path, expected_md5, "complete")
                    return full_path
                logger.warning("File %s exists but MD5 does not match, redownloading", filename)

            # Establish FTP connection
            with ftplib.FTP(host) as ftp:
                ftp.login()  # anonymous login

                # Change to the directory containing the file
                ftp.cwd(os.path.dirname(remote_path))

                # Download the file
                with open(full_path, "wb") as local_file:
                    ftp.retrbinary(f"RETR {os.path.basename(remote_path)}", local_file.write)

            # Verify MD5 checksum if provided
            if expected_md5:
                if self.verify_md5(full_path, expected_md5):
                    logger.info("MD5 checksum verification successful for %s", filename)
                    self._update_download_log(full_path, expected_md5, "complete")
                else:
                    logger.error("MD5 checksum verification failed for %s", filename)
                    self._update_download_log(full_path, expected_md5, "failed_md5")
                    return None

            logger.info("Successfully downloaded: %s", filename)
            return full_path

        except Exception as e:
            logger.error("Failed to download %s: %s", url, e)
            if expected_md5 and "full_path" in locals():
                self._update_download_log(full_path, expected_md5, "failed_download")
            return None

    def download_file(
        self, url: str, filename: str | None = None, subdir: str = "sequencing", expected_md5: str | None = None
    ) -> str | None:
        """
        Download a single file from a given URL.

        Args:
            url: URL of the file to download.
            filename: Optional custom filename.
            subdir: Subdirectory under 'raw' to save the file.
            expected_md5: Expected MD5 checksum for verification.

        Returns:
            Path to the downloaded file or None if download fails.
        """
        try:
            # Check if it's an FTP URL
            if url.startswith("ftp://"):
                return self.download_file_ftp(url, filename, subdir, expected_md5)

            # For HTTP/HTTPS URLs
            if not url.startswith(("http://", "https://")):
                url = f"https://{url}"

            # Extract filename if not provided
            if not filename:
                filename = os.path.basename(urlparse(url).path)

            full_path = os.path.join(self.output_dir, "raw", subdir, filename)

            # Check if file already exists and if we should resume
            downloaded_files = self._get_downloaded_files()
            if self.resume and full_path in downloaded_files:
                file_info = downloaded_files[full_path]
                if file_info["status"] == "complete":
                    if expected_md5 and file_info["md5"] == expected_md5:
                        logger.info("File %s already exists with matching MD5, skipping download", filename)
                        return full_path
                    if os.path.exists(full_path):
                        # Verify the existing file's MD5
                        if expected_md5 and self.verify_md5(full_path, expected_md5):
                            logger.info("File %s already exists with matching MD5, skipping download", filename)
                            return full_path
                        logger.warning("File %s exists but MD5 does not match, redownloading", filename)
                    else:
                        logger.warning("File %s is in log but not found on disk, redownloading", filename)

            # Check if file exists on disk but not in the log
            if self.resume and os.path.exists(full_path) and expected_md5:
                if self.verify_md5(full_path, expected_md5):
                    logger.info("File %s already exists with matching MD5, skipping download", filename)
                    self._update_download_log(full_path, expected_md5, "complete")
                    return full_path
                logger.warning("File %s exists but MD5 does not match, redownloading", filename)

            logger.info("Attempting to download: %s", url)

            # Use requests for HTTP/HTTPS downloads
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(full_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            # Verify MD5 checksum if provided
            if expected_md5:
                if self.verify_md5(full_path, expected_md5):
                    logger.info("MD5 checksum verification successful for %s", filename)
                    self._update_download_log(full_path, expected_md5, "complete")
                else:
                    logger.error("MD5 checksum verification failed for %s", filename)
                    self._update_download_log(full_path, expected_md5, "failed_md5")
                    return None

            logger.info("Successfully downloaded: %s", filename)
            return full_path

        except Exception as e:
            logger.error("Failed to download %s: %s", url, e)
            if expected_md5 and "full_path" in locals():
                self._update_download_log(full_path, expected_md5, "failed_download")
            return None

    def download_dataset(self, accession: str, max_files: int | None = None) -> list[str]:
        """
        Download entire dataset for a given accession.

        Args:
            accession: Study or project accession.
            max_files: Maximum number of files to download (randomly selected).

        Returns:
            List of paths to downloaded files.
        """
        # Search for dataset details
        dataset_info = self.search_dataset(accession)

        # Track downloaded files
        downloaded_files: list[str] = []

        # Prepare a list of all files to download
        files_to_download = []

        # Collect all available files
        for entry in dataset_info:
            if entry.get("fastq_ftp"):
                ftp_urls = entry["fastq_ftp"].split(";")
                # Get MD5 checksums if available
                md5_checksums = []
                if entry.get("fastq_md5"):
                    md5_checksums = entry["fastq_md5"].split(";")

                for i, url in enumerate(ftp_urls):
                    expected_md5 = md5_checksums[i] if i < len(md5_checksums) else None
                    files_to_download.append(
                        {"url": url, "subdir": "sequencing", "filename": None, "md5": expected_md5}
                    )

            if entry.get("submitted_ftp"):
                sample_url = entry["submitted_ftp"].split(";")[0]  # Assume first file is sample metadata
                sample_filename = f"{entry['sample_accession']}.txt"
                files_to_download.append(
                    {"url": sample_url, "subdir": "metadata", "filename": sample_filename, "md5": None}
                )

        # Randomly select files if max_files is specified
        if max_files is not None and max_files < len(files_to_download):
            logger.info("Randomly selecting %d files out of %d available", max_files, len(files_to_download))
            files_to_download = random.sample(files_to_download, max_files)

        # Download selected files
        for file_info in files_to_download:
            downloaded_file = self.download_file(
                url=file_info["url"],
                filename=file_info["filename"],
                subdir=file_info["subdir"],
                expected_md5=file_info["md5"],
            )
            if downloaded_file:
                downloaded_files.append(downloaded_file)

        return downloaded_files

    def verify_md5(self, file_path: str, expected_md5: str) -> bool:
        """
        Verify MD5 checksum of a downloaded file.

        Args:
            file_path: Path to the downloaded file.
            expected_md5: Expected MD5 checksum.

        Returns:
            Boolean indicating if MD5 matches.
        """
        # Calculate MD5 of the file
        md5_hash = hashlib.md5()
        with open(file_path, "rb") as f:
            # Read the file in chunks
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)

        calculated_md5 = md5_hash.hexdigest()

        return calculated_md5 == expected_md5

    def get_download_stats(self) -> dict[str, int]:
        """
        Get statistics about downloaded files.

        Returns:
            Dictionary with statistics.
        """
        downloaded_files = self._get_downloaded_files()

        stats = {"total": len(downloaded_files), "complete": 0, "failed_md5": 0, "failed_download": 0, "other": 0}

        for file_info in downloaded_files.values():
            status = file_info.get("status", "")
            if status == "complete":
                stats["complete"] += 1
            elif status == "failed_md5":
                stats["failed_md5"] += 1
            elif status == "failed_download":
                stats["failed_download"] += 1
            else:
                stats["other"] += 1

        return stats
