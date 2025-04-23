"""
Enhanced ENA Downloader - With comprehensive sample metadata support
"""

import argparse
import csv
import ftplib
import hashlib
import json
import logging
import os
import re
import sys
from urllib.parse import urlparse

import requests

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


def search_ena_dataset(accession):
    """Search for dataset details in ENA."""
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "read_run",
        "query": f"study_accession={accession}",
        "fields": "study_accession,sample_accession,run_accession,fastq_ftp,submitted_ftp,fastq_md5",
        "format": "JSON",
    }

    try:
        logger.info("Searching for dataset: %s", accession)
        response = requests.get(url, params=params)
        response.raise_for_status()

        dataset_info = response.json()
        logger.info("Found %d entries for %s", len(dataset_info), accession)

        return dataset_info
    except requests.RequestException as e:
        logger.error("Error searching dataset: %s", e)
        return []


def get_sample_metadata(sample_accessions):
    """
    Retrieve metadata for a list of sample accessions.

    Args:
        sample_accessions: List of sample accession IDs

    Returns:
        Dictionary mapping sample accessions to their metadata
    """
    if not sample_accessions:
        return {}

    # Get a unique list of sample accessions
    unique_samples = list(set(sample_accessions))
    logger.info("Retrieving metadata for %d unique samples", len(unique_samples))

    # Build a query with all sample accessions (limited batches if needed)
    metadata_by_sample = {}

    # Process in batches of 100 to avoid overly long URLs
    batch_size = 100
    for i in range(0, len(unique_samples), batch_size):
        batch = unique_samples[i : i + batch_size]

        # Build the sample query string
        sample_query = " OR ".join([f'sample_accession="{acc}"' for acc in batch])

        # Build the API request
        url = "https://www.ebi.ac.uk/ena/portal/api/search"
        params = {
            "result": "sample",
            "query": f"({sample_query})",
            "fields": "all",  # Request all available fields
            "format": "JSON",
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()

            # Process each sample's metadata
            for sample_data in response.json():
                sample_acc = sample_data.get("sample_accession")
                if sample_acc:
                    metadata_by_sample[sample_acc] = sample_data

        except requests.RequestException as e:
            logger.error("Error retrieving sample metadata (batch %d): %s", i // batch_size + 1, e)

    logger.info("Retrieved metadata for %d samples", len(metadata_by_sample))
    return metadata_by_sample


def get_file_list(dataset_info, pattern=None, exclude_pattern=None, use_submitted=False):
    """
    Get list of files with priority given to standardized files.

    Args:
        dataset_info: Dataset information from ENA API
        pattern: Regex pattern to include files
        exclude_pattern: Regex pattern to exclude files
        use_submitted: If True, use submitted files instead of standardized fastq

    Returns:
        List of file info dictionaries with duplicates removed
    """
    # Compile regex patterns if provided
    include_regex = re.compile(pattern) if pattern else None
    exclude_regex = re.compile(exclude_pattern) if exclude_pattern else None

    # Dictionary to track files by run_accession to avoid duplicates
    # Key: run_accession
    # Value: list of file information dictionaries
    files_by_run = {}

    for entry in dataset_info:
        run_accession = entry.get("run_accession")
        if not run_accession:
            continue

        # Track files for this run
        if run_accession not in files_by_run:
            files_by_run[run_accession] = []

        # Process files based on preference (standardized by default)
        if not use_submitted and entry.get("fastq_ftp"):
            # Use standardized FASTQ files (default)
            ftp_urls = entry["fastq_ftp"].split(";")
            md5_checksums = entry.get("fastq_md5", "").split(";") if entry.get("fastq_md5") else []

            for i, url_path in enumerate(ftp_urls):
                filename = os.path.basename(urlparse(url_path).path)  # Get original filename
                md5 = md5_checksums[i] if i < len(md5_checksums) else None

                # Ensure URL has proper protocol prefix
                url = f"ftp://{url_path}" if not url_path.startswith(("ftp://", "http://", "https://")) else url_path

                # Apply filters
                if include_regex and not include_regex.search(filename):
                    continue
                if exclude_regex and exclude_regex.search(filename):
                    continue

                files_by_run[run_accession].append(
                    {
                        "filename": filename,
                        "url": url,
                        "type": "fastq",
                        "sample": entry.get("sample_accession"),
                        "run": run_accession,
                        "md5": md5,
                        "size_estimate": "large",
                    }
                )

        elif use_submitted and entry.get("submitted_ftp"):
            # Use original submitted files if requested
            ftp_urls = entry["submitted_ftp"].split(";")

            for url_path in ftp_urls:
                filename = os.path.basename(urlparse(url_path).path)  # Get original filename

                # Ensure URL has proper protocol prefix
                url = f"ftp://{url_path}" if not url_path.startswith(("ftp://", "http://", "https://")) else url_path

                # Apply filters
                if include_regex and not include_regex.search(filename):
                    continue
                if exclude_regex and exclude_regex.search(filename):
                    continue

                # Make a guess about file type based on extension
                file_ext = filename.split(".")[-1].lower() if "." in filename else ""
                is_likely_binary = file_ext in ("bam", "cram", "sra", "gz", "zip", "tar", "fastq", "fq")

                files_by_run[run_accession].append(
                    {
                        "filename": filename,
                        "url": url,
                        "type": "submitted",
                        "sample": entry.get("sample_accession"),
                        "run": run_accession,
                        "md5": None,
                        "size_estimate": "large" if is_likely_binary else "small",
                    }
                )

    # Flatten the dictionary into a list
    all_files = []
    for run_files in files_by_run.values():
        all_files.extend(run_files)

    return all_files


def verify_md5(file_path, expected_md5):
    """
    Verify MD5 checksum of a file.

    Args:
        file_path: Path to the file
        expected_md5: Expected MD5 checksum

    Returns:
        True if MD5 matches, False otherwise
    """
    if not expected_md5:
        logger.warning("No MD5 checksum provided for %s, skipping verification", os.path.basename(file_path))
        return True

    md5_hash = hashlib.md5()
    try:
        with open(file_path, "rb") as f:
            # Read file in chunks to handle large files
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)

        calculated_md5 = md5_hash.hexdigest()

        if calculated_md5 == expected_md5:
            logger.info("MD5 verification successful for %s", os.path.basename(file_path))
            return True
        logger.error("MD5 verification failed for %s", os.path.basename(file_path))
        logger.error("Expected: %s", expected_md5)
        logger.error("Calculated: %s", calculated_md5)
        return False
    except Exception as e:
        logger.error("Error verifying MD5 for %s: %s", file_path, e)
        return False


def check_file_status(file_path, expected_md5):
    """
    Check if a file exists and has the correct MD5.

    Args:
        file_path: Path to the file
        expected_md5: Expected MD5 checksum

    Returns:
        "complete" if file exists and MD5 matches
        "incomplete" if file exists but MD5 doesn't match
        "missing" if file doesn't exist
    """
    if not os.path.exists(file_path):
        return "missing"

    if not expected_md5:
        # If no MD5 is provided, check if file size is > 0
        if os.path.getsize(file_path) > 0:
            return "complete"
        return "incomplete"

    # Verify MD5 if provided
    if verify_md5(file_path, expected_md5):
        return "complete"
    return "incomplete"


def download_file(url, output_dir, expected_md5=None, force=False):
    """
    Download a single file.

    Args:
        url: URL of the file to download
        output_dir: Directory to save the file
        expected_md5: Expected MD5 checksum for verification
        force: If True, download even if file exists

    Returns:
        Path to the downloaded file or None if download fails
    """
    try:
        # Extract filename from URL
        parsed_url = urlparse(url)
        filename = os.path.basename(parsed_url.path)
        output_path = os.path.join(output_dir, filename)

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Check if file already exists and verify it
        if not force and os.path.exists(output_path):
            status = check_file_status(output_path, expected_md5)

            if status == "complete":
                logger.info("File %s already exists and is complete, skipping download", filename)
                return output_path
            if status == "incomplete":
                logger.info("File %s exists but is incomplete or corrupted, redownloading", filename)
            else:
                # Should never reach here as status would be "missing" if not found
                pass

        logger.info("Downloading %s", filename)

        # Download the file
        if url.startswith("ftp://"):
            # FTP download
            ftp_host = parsed_url.netloc
            ftp_path = parsed_url.path

            with ftplib.FTP(ftp_host) as ftp:
                ftp.login()  # Anonymous login
                with open(output_path, "wb") as f:
                    ftp.retrbinary(f"RETR {ftp_path}", f.write)
        else:
            # HTTP(S) download
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(output_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

        # Verify the downloaded file
        if expected_md5 and not verify_md5(output_path, expected_md5):
            logger.error("MD5 verification failed for %s", filename)
            return None

        logger.info("Successfully downloaded %s", filename)
        return output_path

    except Exception as e:
        logger.error("Failed to download %s: %s", url, e)
        if "output_path" in locals() and os.path.exists(output_path):
            # If file was partially downloaded, keep it for potential resuming later
            logger.info("Partial file remains at %s", output_path)
        return None


def save_metadata_to_file(metadata, output_file, format="json"):
    """
    Save sample metadata to a file.

    Args:
        metadata: Dictionary of sample metadata
        output_file: Path to output file
        format: Output format (json or csv)

    Returns:
        True if successful, False otherwise
    """
    try:
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

        if format.lower() == "json":
            with open(output_file, "w") as f:
                json.dump(metadata, f, indent=2)

        elif format.lower() == "csv":
            # Get a union of all fields across all samples
            all_fields = set()
            for sample_data in metadata.values():
                all_fields.update(sample_data.keys())

            field_list = sorted(all_fields)

            with open(output_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=field_list)
                writer.writeheader()
                for sample_data in metadata.values():
                    writer.writerow({k: v for k, v in sample_data.items() if k in field_list})

        logger.info("Saved metadata to %s", output_file)
        return True

    except Exception as e:
        logger.error("Failed to save metadata: %s", e)
        return False


def save_table_to_file(data, headers, output_file):
    """
    Save a table to a text file.

    Args:
        data: List of rows (each row is a list or tuple of values)
        headers: List of column headers
        output_file: Path to output file

    Returns:
        True if successful, False otherwise
    """
    try:
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

        with open(output_file, "w") as f:
            # Calculate column widths based on data
            col_widths = [len(str(h)) for h in headers]
            for row in data:
                for i, val in enumerate(row):
                    col_widths[i] = max(col_widths[i], len(str(val)))

            # Add some padding
            col_widths = [w + 2 for w in col_widths]

            # Write header
            header_row = ""
            for i, header in enumerate(headers):
                header_row += f"{header!s:<{col_widths[i]}}"
            f.write(header_row + "\n")

            # Write separator
            f.write("-" * sum(col_widths) + "\n")

            # Write data rows
            for row in data:
                data_row = ""
                for i, val in enumerate(row):
                    data_row += f"{val!s:<{col_widths[i]}}"
                f.write(data_row + "\n")

        logger.info("Saved table to %s", output_file)
        return True

    except Exception as e:
        logger.error("Failed to save table: %s", e)
        return False


def main():
    """Main function to run the improved ENA downloader."""
    parser = argparse.ArgumentParser(
        description="List, download files, and fetch metadata from the European Nucleotide Archive (ENA)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("accession", help="Study or project accession to search (e.g., PRJNA123456)")
    parser.add_argument(
        "--pattern",
        type=str,
        default=None,
        help="Regex pattern to filter files (only files matching this pattern will be included)",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        default=None,
        help="Regex pattern to exclude files (files matching this pattern will be skipped)",
    )
    parser.add_argument(
        "--use-submitted",
        action="store_true",
        default=False,
        help="Use submitted files instead of standardized FASTQ files",
    )
    parser.add_argument(
        "--download", action="store_true", default=False, help="Download the files (default is to only list them)"
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        default=False,
        help="Force download even if files exist and are complete",
    )
    parser.add_argument("-o", "--output-dir", default="ena_downloads", help="Directory to save downloaded files")
    parser.add_argument("--max-files", type=int, default=None, help="Maximum number of files to download")
    parser.add_argument("--format", choices=["table", "urls"], default="table", help="Output format for listing files")
    parser.add_argument("--metadata", action="store_true", default=False, help="Download and display sample metadata")
    parser.add_argument(
        "--metadata-format", choices=["json", "csv"], default="json", help="Format for saving metadata files"
    )
    parser.add_argument(
        "--metadata-file",
        type=str,
        default=None,
        help="Output file for sample metadata (default: <accession>_metadata.<format> in output dir)",
    )
    parser.add_argument("--save-tables", action="store_true", default=False, help="Save tables to text files")

    args = parser.parse_args()

    # Create project-specific directory structure
    project_dir = os.path.join(args.output_dir, args.accession)
    metadata_dir = os.path.join(project_dir, "metadata")
    data_dir = os.path.join(project_dir, "data")

    # Create directories
    if args.download or args.metadata or args.save_tables:
        os.makedirs(project_dir, exist_ok=True)
        os.makedirs(metadata_dir, exist_ok=True)
        os.makedirs(data_dir, exist_ok=True)

    # Search for dataset details
    dataset_info = search_ena_dataset(args.accession)

    if not dataset_info:
        logger.error("No data found for accession: %s", args.accession)
        return 1

    # Get filtered file list
    files = get_file_list(dataset_info, args.pattern, args.exclude, args.use_submitted)

    # Display results
    if not files:
        logger.info("No files match the specified criteria")
        return 0

    # Apply max files limit if specified
    if args.max_files is not None and args.max_files < len(files):
        logger.info("Limiting to %d files (from %d available)", args.max_files, len(files))
        files = files[: args.max_files]
    else:
        logger.info("Found %d files matching criteria", len(files))

    # Prepare file table data
    file_table_data = []
    for i, file in enumerate(files, 1):
        file_table_data.append((i, file["filename"], file["type"], file["sample"], file["run"]))
    file_headers = ["#", "Filename", "Type", "Sample", "Run"]

    # Display file information
    if args.format == "table":
        # Print table header
        print(f"{'#':<4} {'Filename':<50} {'Type':<10} {'Sample':<15} {'Run':<12}")
        print("-" * 91)

        for i, file in enumerate(files, 1):
            print(f"{i:<4} {file['filename']:<50} {file['type']:<10} {file['sample']:<15} {file['run']:<12}")
    else:
        # Print just the URLs
        for i, file in enumerate(files, 1):
            print(f"{i}. {file['url']}")

    # Save file table if requested
    if args.save_tables:
        file_table_path = os.path.join(metadata_dir, f"{args.accession}_files.txt")
        save_table_to_file(file_table_data, file_headers, file_table_path)
        print(f"\nFile list saved to: {file_table_path}")

    # Download files if requested
    if args.download:
        logger.info("Downloading %d files to %s", len(files), data_dir)

        # Pre-check files to see which need downloading
        if not args.force_download:
            status_summary = {"complete": 0, "incomplete": 0, "missing": 0}
            for file in files:
                file_path = os.path.join(data_dir, file["filename"])
                status = check_file_status(file_path, file.get("md5"))
                status_summary[status] += 1

            if status_summary["complete"] > 0:
                logger.info("%d of %d files already complete", status_summary["complete"], len(files))
            if status_summary["incomplete"] > 0:
                logger.info("%d of %d files exist but need redownloading", status_summary["incomplete"], len(files))
            if status_summary["missing"] > 0:
                logger.info("%d of %d files need to be downloaded", status_summary["missing"], len(files))

            # If all files are already complete, we can skip download
            if status_summary["complete"] == len(files):
                logger.info("All files already exist and are complete, skipping download")
                # Continue with metadata processing

        # Download files
        downloaded_files = []
        for file in files:
            downloaded_file = download_file(file["url"], data_dir, file.get("md5"), force=args.force_download)
            if downloaded_file:
                downloaded_files.append(downloaded_file)

        logger.info("Successfully downloaded %d of %d files", len(downloaded_files), len(files))

    # Process metadata if requested
    if args.metadata:
        # Extract sample accessions from files
        sample_accessions = [file["sample"] for file in files if file.get("sample")]

        # Get metadata for all samples
        metadata = get_sample_metadata(sample_accessions)

        if metadata:
            # Determine metadata output file
            if args.metadata_file:
                metadata_file = os.path.join(metadata_dir, args.metadata_file)
            else:
                extension = ".json" if args.metadata_format == "json" else ".csv"
                metadata_file = os.path.join(metadata_dir, f"{args.accession}_metadata{extension}")

            # Save metadata to file
            save_metadata_to_file(metadata, metadata_file, args.metadata_format)

            # Prepare metadata table data
            metadata_table_data = []
            for i, (sample_id, sample_data) in enumerate(sorted(metadata.items()), 1):
                organism = sample_data.get("scientific_name", "Unknown")
                description = sample_data.get("description", "No description")

                # Truncate long descriptions for display but not for saving
                display_description = description
                if len(display_description) > 40:
                    display_description = display_description[:37] + "..."

                metadata_table_data.append((i, sample_id, organism, description))

            metadata_headers = ["#", "Sample ID", "Organism", "Description"]

            # Display summary of metadata
            print("\nSample Metadata Summary:")
            print(f"{'#':<4} {'Sample ID':<15} {'Organism':<25} {'Description':<40}")
            print("-" * 84)

            for i, (sample_id, sample_data) in enumerate(sorted(metadata.items()), 1):
                organism = sample_data.get("scientific_name", "Unknown")
                description = sample_data.get("description", "No description")

                # Truncate long descriptions
                if len(description) > 40:
                    description = description[:37] + "..."

                print(f"{i:<4} {sample_id:<15} {organism:<25} {description:<40}")

            print(f"\nFull metadata saved to: {metadata_file}")

            # Save metadata table if requested
            if args.save_tables:
                metadata_table_path = os.path.join(metadata_dir, f"{args.accession}_samples.txt")
                save_table_to_file(metadata_table_data, metadata_headers, metadata_table_path)
                print(f"Sample summary saved to: {metadata_table_path}")
        else:
            logger.warning("No metadata retrieved for samples")

    return 0


if __name__ == "__main__":
    sys.exit(main())
