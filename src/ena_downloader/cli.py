#!/usr/bin/env python3
"""
Enhanced ENA Downloader - With comprehensive accession support and direct assembly access
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
import time
from urllib.parse import urlparse

import requests

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)

# Define standard directory structure for genomics data based on accession type
DIR_STRUCTURE = {
    # Common directories for all accession types
    "raw_data": "raw-data",  # Raw sequencing data (FASTQ)
    "assemblies": "assemblies",  # Genome assemblies (FASTA)
    "analyses": "analyses",  # Analysis results
    "metadata": "metadata",  # Metadata files
    "annotations": "annotations",  # Genome annotations (GFF, GTF)
    "alignments": "alignments",  # Alignment files (BAM, CRAM)
    "variants": "variants",  # Variant calls (VCF)
    "runs": "runs",  # Run information
    "experiments": "experiments",  # Experiment information
    "samples": "samples",  # Sample-specific information
    "sequences": "sequences",  # Assembled/Annotated Sequences
    "proteins": "proteins",  # Protein Coding Sequences
}

# Define accession patterns
ACCESSION_PATTERNS = {
    "project": [
        r"^PRJ[EDN][A-Z][0-9]+$",  # PRJEB12345, PRJNA12345, etc.
        r"^[EDS]RP[0-9]{6,}$",  # ERP123456, SRP123456, etc.
    ],
    "biosample": [
        r"^SAM[EDN][A-Z]?[0-9]+$",  # SAMEA123456, SAMN123456, etc.
        r"^[EDS]RS[0-9]{6,}$",  # ERS123456, SRS123456, etc.
    ],
    "experiment": [
        r"^[EDS]RX[0-9]{6,}$"  # ERX123456, SRX123456, etc.
    ],
    "run": [
        r"^[EDS]RR[0-9]{6,}$"  # ERR123456, SRR123456, etc.
    ],
    "analysis": [
        r"^[EDS]RZ[0-9]{6,}$"  # ERZ123456, SRZ123456, etc.
    ],
    "assembly": [
        r"^GC[AF]_[0-9]{9}\.[0-9]+$"  # GCA_123456789.1, GCF_123456789.1, etc.
    ],
    "sequence": [
        r"^[A-Z]{1}[0-9]{5}\.[0-9]+$",  # A12345.1, etc.
        r"^[A-Z]{2}[0-9]{6}\.[0-9]+$",  # AB123456.1, etc.
        r"^[A-Z]{2}[0-9]{8}$",  # AB12345678, etc.
        r"^[A-Z]{4}[0-9]{2}S?[0-9]{6,8}$",  # ABCD01123456, etc.
        r"^[A-Z]{6}[0-9]{2}S?[0-9]{7,9}$",  # ABCDEF011234567, etc.
    ],
    "protein": [
        r"^[A-Z]{3}[0-9]{5}\.[0-9]+$",  # ABC12345.1, etc.
        r"^[A-Z]{3}[0-9]{7}\.[0-9]+$",  # ABC1234567.1, etc.
    ],
}


def identify_accession_type(accession):
    """
    Identify the type of accession based on its format.

    Args:
        accession: Accession to check

    Returns:
        Tuple of (accession_type, True/False) where accession_type is one of:
        'project', 'biosample', 'experiment', 'run', 'analysis', 'assembly',
        'sequence', 'protein', or None if unrecognized
    """
    for acc_type, patterns in ACCESSION_PATTERNS.items():
        for pattern in patterns:
            if re.match(pattern, accession):
                return acc_type, True

    return None, False


def is_assembly_accession(accession):
    """
    Check if the accession is an assembly accession.

    Args:
        accession: Accession to check

    Returns:
        True if assembly accession, False otherwise
    """
    acc_type, is_valid = identify_accession_type(accession)
    return acc_type == "assembly" and is_valid


def is_study_accession(accession):
    """
    Check if the accession is a study/project accession.

    Args:
        accession: Accession to check

    Returns:
        True if study accession, False otherwise
    """
    acc_type, is_valid = identify_accession_type(accession)
    return acc_type == "project" and is_valid


def search_ena_dataset(accession):
    """Search for dataset details in ENA."""
    url = "https://www.ebi.ac.uk/ena/portal/api/search"

    # Modify the query based on accession type
    acc_type, _ = identify_accession_type(accession)

    if acc_type == "project":
        query_field = "study_accession"
    elif acc_type == "biosample":
        query_field = "sample_accession"
    elif acc_type == "experiment":
        query_field = "experiment_accession"
    elif acc_type == "run":
        query_field = "run_accession"
    elif acc_type == "analysis":
        query_field = "analysis_accession"
    else:
        # Default to study_accession for backward compatibility
        query_field = "study_accession"

    params = {
        "result": "read_run",
        "query": f"{query_field}={accession}",
        "fields": "study_accession,sample_accession,experiment_accession,run_accession,fastq_ftp,submitted_ftp,fastq_md5",
        "format": "JSON",
    }

    try:
        logger.info("Searching for dataset: %s (type: %s)", accession, acc_type)
        response = requests.get(url, params=params)
        response.raise_for_status()

        dataset_info = response.json()
        logger.info("Found %d entries for %s", len(dataset_info), accession)

        return dataset_info
    except requests.RequestException as e:
        logger.error("Error searching dataset: %s", e)
        return []


def search_assembly(accession):
    """
    Simple assembly search that creates a minimal info object for browser API access.

    Args:
        accession: Assembly accession (e.g., GCA_031216745.1)

    Returns:
        Minimal assembly info object sufficient for browser API download
    """
    logger.info(f"Creating minimal assembly object for direct download: {accession}")

    # Create a minimal assembly info object with just the accession
    return [{"assembly_accession": accession, "assembly_name": accession, "assembly_title": f"Assembly {accession}"}]


def search_sequence(accession):
    """
    Search for sequence details in ENA.

    Args:
        accession: Sequence accession

    Returns:
        Sequence information from ENA API
    """
    # Remove version from accession for broader search if it has one
    base_accession = accession.split(".")[0] if "." in accession else accession

    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "sequence",
        "query": f"sequence_accession={base_accession}",
        "fields": "sequence_accession,sequence_version,sequence_name,sequence_length,study_accession,sample_accession,tax_id,scientific_name",
        "format": "JSON",
    }

    try:
        logger.info("Searching for sequence: %s", accession)
        response = requests.get(url, params=params)
        response.raise_for_status()

        sequences = response.json()

        # Filter to match the specific version if provided
        if "." in accession:
            sequences = [
                s
                for s in sequences
                if s.get("sequence_accession") == base_accession
                and str(s.get("sequence_version")) == accession.split(".")[1]
            ]

        if sequences:
            logger.info("Found sequence: %s", accession)
            return sequences

        logger.error("No sequence found for accession: %s", accession)
        return []

    except requests.RequestException as e:
        logger.error("Error searching sequence: %s", e)
        return []


def search_protein(accession):
    """
    Search for protein sequence details in ENA.

    Args:
        accession: Protein accession

    Returns:
        Protein sequence information from ENA API
    """
    # Remove version from accession for broader search if it has one
    base_accession = accession.split(".")[0] if "." in accession else accession

    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "coding",
        "query": f"coding_accession={base_accession}",
        "fields": "coding_accession,coding_version,protein_id,coding_product,study_accession,sample_accession,tax_id,scientific_name",
        "format": "JSON",
    }

    try:
        logger.info("Searching for protein: %s", accession)
        response = requests.get(url, params=params)
        response.raise_for_status()

        proteins = response.json()

        # Filter to match the specific version if provided
        if "." in accession:
            proteins = [
                p
                for p in proteins
                if p.get("coding_accession") == base_accession
                and str(p.get("coding_version")) == accession.split(".")[1]
            ]

        if proteins:
            logger.info("Found protein: %s", proteins[0].get("coding_product", "Unknown protein"))
            return proteins

        logger.error("No protein found for accession: %s", accession)
        return []

    except requests.RequestException as e:
        logger.error("Error searching protein: %s", e)
        return []


def get_assembly_files(assembly_info, include_embl=True):
    """
    Get list of files for an assembly using direct browser API links.

    Args:
        assembly_info: Assembly information (can be minimal with just accession)

    Returns:
        List of file info dictionaries
    """
    files = []

    for entry in assembly_info:
        assembly_accession = entry.get("assembly_accession")
        if not assembly_accession:
            continue

        # Use direct browser API links for reliable access
        # Main genomic FASTA file
        fasta_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{assembly_accession}?download=true&gzip=true"
        files.append(
            {
                "filename": f"{assembly_accession}_genomic.fasta.gz",
                "url": fasta_url,
                "type": "genomic_fasta",
                "description": "Complete Genomic FASTA",
                "assembly": assembly_accession,
                "sample": entry.get("sample_accession"),
                "project": entry.get("study_accession"),
                "md5": None,
                "size_estimate": "large",
            }
        )

    if include_embl:
        embl_url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{assembly_accession}?download=true&gzip=true"
        files.append(
            {
                "filename": f"{assembly_accession}_assembly.embl.gz",
                "url": embl_url,
                "type": "assembly_annotation",  # Changed from "assembly_embl"
                "description": "Assembly in EMBL format",
                "assembly": assembly_accession,
                "sample": entry.get("sample_accession"),
                "project": entry.get("study_accession"),
                "md5": None,
                "size_estimate": "large",
            }
        )

    return files


def get_sequence_files(sequence_info):
    """
    Get list of files for a sequence.

    Args:
        sequence_info: Sequence information from ENA API

    Returns:
        List of file info dictionaries
    """
    files = []

    for entry in sequence_info:
        sequence_accession = entry.get("sequence_accession")
        if not sequence_accession:
            continue

        sequence_version = entry.get("sequence_version", "1")
        full_accession = f"{sequence_accession}.{sequence_version}"

        # Construct ENA file URLs for the sequence
        # Typically, sequence data can be downloaded in FASTA and EMBL format

        # EMBL format (flat file)
        embl_url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{full_accession}"
        files.append(
            {
                "filename": f"{full_accession}.embl",
                "url": embl_url,
                "type": "sequence_embl",
                "description": "Sequence in EMBL format",
                "sequence": full_accession,
                "sample": entry.get("sample_accession"),
                "project": entry.get("study_accession"),
                "md5": None,
                "size_estimate": "medium",
            }
        )

        # FASTA format
        fasta_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{full_accession}"
        files.append(
            {
                "filename": f"{full_accession}.fasta",
                "url": fasta_url,
                "type": "sequence_fasta",
                "description": "Sequence in FASTA format",
                "sequence": full_accession,
                "sample": entry.get("sample_accession"),
                "project": entry.get("study_accession"),
                "md5": None,
                "size_estimate": "medium",
            }
        )

    return files


def get_protein_files(protein_info):
    """
    Get list of files for a protein.

    Args:
        protein_info: Protein information from ENA API

    Returns:
        List of file info dictionaries
    """
    files = []

    for entry in protein_info:
        coding_accession = entry.get("coding_accession")
        if not coding_accession:
            continue

        coding_version = entry.get("coding_version", "1")
        full_accession = f"{coding_accession}.{coding_version}"

        # Construct ENA file URLs for the protein sequence

        # FASTA format
        fasta_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{full_accession}"
        files.append(
            {
                "filename": f"{full_accession}.fasta",
                "url": fasta_url,
                "type": "protein_fasta",
                "description": "Protein sequence in FASTA format",
                "protein": full_accession,
                "sample": entry.get("sample_accession"),
                "project": entry.get("study_accession"),
                "md5": None,
                "size_estimate": "small",
            }
        )

    return files


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


def get_project_metadata(accession):
    """
    Retrieve metadata for a project accession.

    Args:
        accession: Project accession ID

    Returns:
        Dictionary with project metadata
    """
    if not accession:
        return {}

    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "study",
        "query": f'study_accession="{accession}"',
        "fields": "all",
        "format": "JSON",
    }

    try:
        logger.info("Retrieving metadata for project: %s", accession)
        response = requests.get(url, params=params)
        response.raise_for_status()

        project_data = response.json()
        if project_data:
            return project_data[0]

        logger.warning("No metadata found for project: %s", accession)
        return {}

    except requests.RequestException as e:
        logger.error("Error retrieving project metadata: %s", e)
        return {}


def get_run_metadata(accession):
    """
    Retrieve metadata for a run accession.

    Args:
        accession: Run accession ID

    Returns:
        Dictionary with run metadata
    """
    if not accession:
        return {}

    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "read_run",
        "query": f'run_accession="{accession}"',
        "fields": "all",
        "format": "JSON",
    }

    try:
        logger.info("Retrieving metadata for run: %s", accession)
        response = requests.get(url, params=params)
        response.raise_for_status()

        run_data = response.json()
        if run_data:
            return run_data[0]

        logger.warning("No metadata found for run: %s", accession)
        return {}

    except requests.RequestException as e:
        logger.error("Error retrieving run metadata: %s", e)
        return {}


def get_experiment_metadata(accession):
    """
    Retrieve metadata for an experiment accession.

    Args:
        accession: Experiment accession ID

    Returns:
        Dictionary with experiment metadata
    """
    if not accession:
        return {}

    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "read_experiment",
        "query": f'experiment_accession="{accession}"',
        "fields": "all",
        "format": "JSON",
    }

    try:
        logger.info("Retrieving metadata for experiment: %s", accession)
        response = requests.get(url, params=params)
        response.raise_for_status()

        experiment_data = response.json()
        if experiment_data:
            return experiment_data[0]

        logger.warning("No metadata found for experiment: %s", accession)
        return {}

    except requests.RequestException as e:
        logger.error("Error retrieving experiment metadata: %s", e)
        return {}


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
                        "experiment": entry.get("experiment_accession"),
                        "run": run_accession,
                        "project": entry.get("study_accession"),
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

                # Determine file type based on extension
                file_type = "submitted"
                if file_ext in ("bam", "cram"):
                    file_type = "alignment"
                elif file_ext in ("vcf", "vcf.gz"):
                    file_type = "variant"
                elif file_ext in ("fastq", "fq", "fastq.gz", "fq.gz"):
                    file_type = "fastq"
                elif file_ext in ("sra"):
                    file_type = "sra"

                files_by_run[run_accession].append(
                    {
                        "filename": filename,
                        "url": url,
                        "type": file_type,
                        "sample": entry.get("sample_accession"),
                        "experiment": entry.get("experiment_accession"),
                        "run": run_accession,
                        "project": entry.get("study_accession"),
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

        # If URL contains query parameters, they won't be in the filename
        if not filename and "?" in url:
            # Extract the base filename before the query parameters
            path = parsed_url.path
            if path:
                filename = os.path.basename(path)
            else:
                # Generate a generic filename using the accession from the URL if possible
                parts = url.split("/")
                for part in parts:
                    if part.startswith("GCA_") or part.startswith("GCF_"):
                        filename = f"{part}.data"
                        break
                if not filename:
                    # Last resort: use a timestamp
                    filename = f"download_{int(time.time())}"

        # Add file extension for browser API downloads if needed
        if "download=true" in url and "gzip=true" in url:
            if "fasta" in url and not filename.endswith(".gz"):
                filename = f"{filename}.fasta.gz"
            elif "embl" in url and not filename.endswith(".gz"):
                filename = f"{filename}.embl.gz"

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


def get_output_directory(base_dir, accession, file_type, accession_type):
    """
    Get the appropriate output directory for a file type based on accession type.

    Args:
        base_dir: Base output directory
        accession: Accession of the project/study/sample/etc.
        file_type: Type of file
        accession_type: Type of accession (project, run, assembly, etc.)

    Returns:
        Path to the appropriate output directory
    """
    # Create project directory - base it on the accession
    project_dir = os.path.join(base_dir, accession)

    # For specific accession types, organize differently
    if accession_type == "project":
        # Project accessions get standard organization
        if file_type in ["fastq", "fastq.gz", "fq", "fq.gz", "sra"]:
            return os.path.join(project_dir, DIR_STRUCTURE["raw_data"])
        if file_type in ["genomic_fasta", "cds_fasta", "protein_fasta", "rna_fasta"]:
            return os.path.join(project_dir, DIR_STRUCTURE["assemblies"])
        if file_type in ["genomic_gff", "genomic_genbank"]:
            return os.path.join(project_dir, DIR_STRUCTURE["annotations"])
        if file_type in ["assembly_report", "assembly_stats"]:
            return os.path.join(project_dir, DIR_STRUCTURE["metadata"], "assembly_info")
        if file_type in ["alignment", "bam", "cram"]:
            return os.path.join(project_dir, DIR_STRUCTURE["alignments"])
        if file_type in ["variant", "vcf", "vcf.gz"]:
            return os.path.join(project_dir, DIR_STRUCTURE["variants"])
        # Default directory for unknown file types
        return os.path.join(project_dir, "other_files")

    if accession_type == "run":
        # Run accessions focus on raw data
        if file_type in ["fastq", "fastq.gz", "fq", "fq.gz", "sra"]:
            return os.path.join(project_dir, DIR_STRUCTURE["raw_data"])
        if file_type in ["alignment", "bam", "cram"]:
            return os.path.join(project_dir, DIR_STRUCTURE["alignments"])
        return os.path.join(project_dir, DIR_STRUCTURE["runs"])

    if accession_type == "experiment":
        # Experiment accessions
        if file_type in ["fastq", "fastq.gz", "fq", "fq.gz", "sra"]:
            return os.path.join(project_dir, DIR_STRUCTURE["raw_data"])
        return os.path.join(project_dir, DIR_STRUCTURE["experiments"])

    if accession_type == "assembly":
        # Assembly accessions focus on assembled data
        if file_type in ["genomic_fasta", "cds_fasta", "protein_fasta", "rna_fasta"]:
            return os.path.join(project_dir, DIR_STRUCTURE["assemblies"])
        if file_type in ["genomic_gff", "genomic_genbank", "assembly_annotation"]:  # Add assembly_annotation here
            return os.path.join(project_dir, DIR_STRUCTURE["annotations"])
        if file_type in ["assembly_report", "assembly_stats"]:  # Remove assembly_embl from here
            return os.path.join(project_dir, DIR_STRUCTURE["metadata"], "assembly_info")
        return os.path.join(project_dir, DIR_STRUCTURE["assemblies"], "other_files")

    if accession_type == "biosample":
        # Biosample accessions focus on sample data
        return os.path.join(project_dir, DIR_STRUCTURE["samples"])

    if accession_type == "sequence":
        # Sequence accessions
        return os.path.join(project_dir, DIR_STRUCTURE["sequences"])

    if accession_type == "protein":
        # Protein accessions
        return os.path.join(project_dir, DIR_STRUCTURE["proteins"])

    if accession_type == "analysis":
        # Analysis accessions
        return os.path.join(project_dir, DIR_STRUCTURE["analyses"])

    # Unknown accession type, use standard mapping
    if file_type in ["fastq", "fastq.gz", "fq", "fq.gz", "sra"]:
        return os.path.join(project_dir, DIR_STRUCTURE["raw_data"])
    if file_type in ["genomic_fasta", "cds_fasta", "protein_fasta", "rna_fasta"]:
        return os.path.join(project_dir, DIR_STRUCTURE["assemblies"])
    if file_type in ["genomic_gff", "genomic_genbank"]:
        return os.path.join(project_dir, DIR_STRUCTURE["annotations"])
    if file_type in ["assembly_report", "assembly_stats"]:
        return os.path.join(project_dir, DIR_STRUCTURE["metadata"], "assembly_info")
    if file_type in ["alignment", "bam", "cram"]:
        return os.path.join(project_dir, DIR_STRUCTURE["alignments"])
    if file_type in ["variant", "vcf", "vcf.gz"]:
        return os.path.join(project_dir, DIR_STRUCTURE["variants"])
    # Default directory for unknown file types
    return os.path.join(project_dir, "other_files")


def determine_project_accession(info_dict):
    """
    Determine the project accession from various sources.

    Args:
        info_dict: Dictionary that might contain project information

    Returns:
        Project accession or None if not found
    """
    # Try different possible keys for project accession
    for key in ["project", "study_accession", "project_accession"]:
        if info_dict.get(key):
            return info_dict[key]

    return None


def main():
    """Main function to run the enhanced ENA downloader with comprehensive accession support."""
    parser = argparse.ArgumentParser(
        description="Enhanced tool to download files and fetch metadata from the European Nucleotide Archive (ENA) with support for all accession types",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("accession", help="Accession to search (e.g., PRJNA123456, GCA_002870075.4, ERR123456, etc.)")
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
        help="Use submitted files instead of standardized FASTQ files (for study/project accessions)",
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
        help="Output file for metadata (default: <accession>_metadata.<format> in metadata dir)",
    )
    parser.add_argument("--save-tables", action="store_true", default=False, help="Save tables to text files")
    parser.add_argument(
        "--file-types",
        type=str,
        default=None,
        help="Comma-separated list of file types to download (e.g., genomic_fasta,genomic_gff)",
    )
    parser.add_argument("--no-embl", action="store_true", default=False, help="Skip downloading EMBL format files")

    args = parser.parse_args()

    # Identify accession type
    accession_type, valid_accession = identify_accession_type(args.accession)

    if not valid_accession:
        logger.error("Unrecognized accession format: %s", args.accession)
        logger.error("Please provide a valid ENA accession number")
        return 1

    logger.info("Processing %s accession: %s", accession_type, args.accession)

    # Create base directory
    base_dir = args.output_dir
    os.makedirs(base_dir, exist_ok=True)

    # Initialize variables
    files = []
    metadata = {}
    project_accession = None

    # Process based on accession type
    if accession_type == "project":
        # Process project/study accession
        dataset_info = search_ena_dataset(args.accession)

        if not dataset_info:
            logger.error("No data found for study accession: %s", args.accession)
            return 1

        # Get file list
        files = get_file_list(dataset_info, args.pattern, args.exclude, args.use_submitted)
        project_accession = args.accession

    elif accession_type == "assembly":
        # Process assembly accession - use simplified direct approach
        # Skip the ENA API search to avoid 400 errors for certain assemblies
        assembly_info = search_assembly(args.accession)

        # Get assembly files using browser API links
        files = get_assembly_files(assembly_info)

        # Filter by file types if specified
        if args.file_types:
            file_types = [t.strip() for t in args.file_types.split(",")]
            files = [f for f in files if f["type"] in file_types]
            logger.info("Filtered to %d files based on specified file types", len(files))

        # Use the assembly accession as the project accession
        project_accession = args.accession.split(".")[0] if "." in args.accession else args.accession

    elif accession_type == "run":
        # Process run accession
        run_info = search_ena_dataset(args.accession)

        if not run_info:
            logger.error("No data found for run accession: %s", args.accession)
            return 1

        # Get file list
        files = get_file_list(run_info, args.pattern, args.exclude, args.use_submitted)

        # Try to get project accession from run info
        if run_info and run_info[0].get("study_accession"):
            project_accession = run_info[0].get("study_accession")
        else:
            project_accession = args.accession

    elif accession_type == "experiment":
        # Process experiment accession
        experiment_info = search_ena_dataset(args.accession)

        if not experiment_info:
            logger.error("No data found for experiment accession: %s", args.accession)
            return 1

        # Get file list
        files = get_file_list(experiment_info, args.pattern, args.exclude, args.use_submitted)

        # Try to get project accession from experiment info
        if experiment_info and experiment_info[0].get("study_accession"):
            project_accession = experiment_info[0].get("study_accession")
        else:
            project_accession = args.accession

    elif accession_type == "biosample":
        # Process biosample/sample accession
        # For biosample accessions, we'll fetch metadata only
        sample_info = get_sample_metadata([args.accession])

        if not sample_info:
            logger.error("No data found for sample accession: %s", args.accession)
            return 1

        metadata = sample_info

        # Try to find runs associated with this sample
        sample_runs = search_ena_dataset(args.accession)
        if sample_runs:
            # Get file list from runs associated with this sample
            files = get_file_list(sample_runs, args.pattern, args.exclude, args.use_submitted)
            logger.info("Found %d files from runs associated with sample %s", len(files), args.accession)

            # Try to get project accession from sample runs
            if sample_runs and sample_runs[0].get("study_accession"):
                project_accession = sample_runs[0].get("study_accession")

        # If no project accession found, use the sample accession
        if not project_accession:
            project_accession = args.accession

    elif accession_type == "sequence":
        # Process sequence accession
        sequence_info = search_sequence(args.accession)

        if not sequence_info:
            logger.error("No sequence found for accession: %s", args.accession)
            return 1

        # Get sequence files
        files = get_sequence_files(sequence_info)

        # Use sequence info as metadata
        metadata = sequence_info[0] if sequence_info else {}

        # Try to get project accession from sequence info
        if sequence_info and sequence_info[0].get("study_accession"):
            project_accession = sequence_info[0].get("study_accession")
        else:
            project_accession = args.accession.split(".")[0] if "." in args.accession else args.accession

    elif accession_type == "protein":
        # Process protein accession
        protein_info = search_protein(args.accession)

        if not protein_info:
            logger.error("No protein found for accession: %s", args.accession)
            return 1

        # Get protein files
        files = get_protein_files(protein_info)

        # Use protein info as metadata
        metadata = protein_info[0] if protein_info else {}

        # Try to get project accession from protein info
        if protein_info and protein_info[0].get("study_accession"):
            project_accession = protein_info[0].get("study_accession")
        else:
            project_accession = args.accession.split(".")[0] if "." in args.accession else args.accession

    else:
        # For other accession types, default to simple search and listing
        logger.warning("Detailed support for %s accession type is limited", accession_type)
        dataset_info = search_ena_dataset(args.accession)

        if dataset_info:
            files = get_file_list(dataset_info, args.pattern, args.exclude, args.use_submitted)

            # Try to get project accession
            if dataset_info and dataset_info[0].get("study_accession"):
                project_accession = dataset_info[0].get("study_accession")

        # If no project accession found, use the given accession
        if not project_accession:
            project_accession = args.accession

    # Display results
    if not files and not metadata:
        logger.info("No files or metadata found matching the specified criteria")
        return 0

    # Apply max files limit if specified and files were found
    if files:
        if args.max_files is not None and args.max_files < len(files):
            logger.info("Limiting to %d files (from %d available)", args.max_files, len(files))
            files = files[: args.max_files]
        else:
            logger.info("Found %d files matching criteria", len(files))

        # Prepare file table data
        file_table_data = []

        # Adjust header based on accession type
        if accession_type == "assembly":
            file_headers = ["#", "Filename", "Type", "Description", "Assembly"]
            for i, file in enumerate(files, 1):
                file_table_data.append(
                    (i, file["filename"], file["type"], file.get("description", ""), file.get("assembly", ""))
                )
        elif accession_type in ["sequence", "protein"]:
            file_headers = ["#", "Filename", "Type", "Description", "Accession"]
            for i, file in enumerate(files, 1):
                file_table_data.append(
                    (i, file["filename"], file["type"], file.get("description", ""), file.get(accession_type, ""))
                )
        else:
            # Default format for project, run, experiment, etc.
            file_headers = ["#", "Filename", "Type", "Sample", "Run"]
            for i, file in enumerate(files, 1):
                file_table_data.append((i, file["filename"], file["type"], file.get("sample", ""), file.get("run", "")))

        # Display file information
        if args.format == "table":
            # Print table header
            header_row = " ".join([f"{h:<{15 if i > 0 else 4}}" for i, h in enumerate(file_headers)])
            print(header_row)
            print("-" * (15 * (len(file_headers) - 1) + 4))

            for row in file_table_data:
                data_row = " ".join([f"{val!s:<{15 if i > 0 else 4}}" for i, val in enumerate(row)])
                print(data_row)
        else:
            # Print just the URLs
            for i, file in enumerate(files, 1):
                print(f"{i}. {file['url']}")

    # Create standard directory structure
    project_dir = os.path.join(base_dir, project_accession)
    # project_dir = base_dir
    metadata_dir = os.path.join(project_dir, DIR_STRUCTURE["metadata"])

    # Create directories if needed for downloads or metadata
    if args.download or args.metadata or args.save_tables:
        os.makedirs(project_dir, exist_ok=True)
        if args.metadata or args.save_tables:
            os.makedirs(metadata_dir, exist_ok=True)

    # Save file table if requested
    if args.save_tables and files:
        file_table_path = os.path.join(metadata_dir, f"{args.accession}_files.txt")
        save_table_to_file(file_table_data, file_headers, file_table_path)
        print(f"\nFile list saved to: {file_table_path}")

    # Download files if requested
    if args.download and files:
        logger.info("Downloading %d files", len(files))

        # Download files
        downloaded_files = []
        for file in files:
            # Determine appropriate output directory based on file type and accession type
            output_dir = get_output_directory(base_dir, project_accession, file["type"], accession_type)

            # Download the file
            downloaded_file = download_file(file["url"], output_dir, file.get("md5"), force=args.force_download)
            if downloaded_file:
                downloaded_files.append(downloaded_file)
                logger.info("Saved to: %s", downloaded_file)

        logger.info("Successfully downloaded %d of %d files", len(downloaded_files), len(files))

    # Process metadata if requested or if it's a metadata-only accession type
    if args.metadata or accession_type == "biosample":
        # Handle metadata based on accession type
        if accession_type == "project":
            # Get project metadata and samples
            sample_accessions = [file["sample"] for file in files if file.get("sample")]
            sample_metadata = get_sample_metadata(sample_accessions)
            project_metadata = get_project_metadata(args.accession)

            metadata = {"project": project_metadata, "samples": sample_metadata}
        elif accession_type == "run":
            # Get run metadata
            run_metadata = get_run_metadata(args.accession)

            if run_metadata:
                # If run has a sample, get sample metadata too
                sample_acc = run_metadata.get("sample_accession")
                sample_metadata = get_sample_metadata([sample_acc]) if sample_acc else {}

                metadata = {
                    "run": run_metadata,
                    "sample": sample_metadata.get(sample_acc) if sample_acc in sample_metadata else None,
                }
        elif accession_type == "experiment":
            # Get experiment metadata
            experiment_metadata = get_experiment_metadata(args.accession)

            if experiment_metadata:
                # If experiment has a sample, get sample metadata too
                sample_acc = experiment_metadata.get("sample_accession")
                sample_metadata = get_sample_metadata([sample_acc]) if sample_acc else {}

                metadata = {
                    "experiment": experiment_metadata,
                    "sample": sample_metadata.get(sample_acc) if sample_acc in sample_metadata else None,
                }
        elif accession_type == "assembly":
            # Use assembly information as metadata
            metadata = assembly_info[0] if assembly_info else {}

        # For biosample, sequence, protein - metadata was already fetched above

        if metadata:
            # Determine metadata output file
            if args.metadata_file:
                metadata_file = os.path.join(metadata_dir, args.metadata_file)
            else:
                extension = ".json" if args.metadata_format == "json" else ".csv"
                metadata_file = os.path.join(metadata_dir, f"{args.accession}_metadata{extension}")

            # Save metadata to file
            save_metadata_to_file(metadata, metadata_file, args.metadata_format)

            print(f"\nMetadata saved to: {metadata_file}")

            # Display metadata summary based on accession type
            if accession_type == "project" and "samples" in metadata and metadata["samples"]:
                # Display summary of sample metadata for project accession
                print("\nSample Metadata Summary:")
                headers = ["#", "Sample ID", "Organism", "Description"]
                print(f"{headers[0]:<4} {headers[1]:<15} {headers[2]:<25} {headers[3]:<40}")
                print("-" * 84)

                sample_table_data = []
                for i, (sample_id, sample_data) in enumerate(sorted(metadata["samples"].items()), 1):
                    organism = sample_data.get("scientific_name", "Unknown")
                    description = sample_data.get("description", "No description")

                    # Truncate long descriptions
                    display_description = description
                    if len(display_description) > 40:
                        display_description = display_description[:37] + "..."

                    sample_table_data.append((i, sample_id, organism, description))
                    print(f"{i:<4} {sample_id:<15} {organism:<25} {display_description:<40}")

                # Save sample table if requested
                if args.save_tables:
                    sample_table_path = os.path.join(metadata_dir, f"{args.accession}_samples.txt")
                    save_table_to_file(sample_table_data, headers, sample_table_path)
                    print(f"Sample summary saved to: {sample_table_path}")

            elif accession_type == "biosample":
                # Display biosample metadata
                print("\nBiosample Metadata Summary:")
                if isinstance(metadata, dict) and args.accession in metadata:
                    sample_data = metadata[args.accession]
                    for key, value in sorted(sample_data.items()):
                        if value and key not in ["sample_xml"]:  # Skip XML and empty values
                            print(f"{key}: {value}")
                else:
                    # Handle flat metadata structure
                    for key, value in sorted(metadata.items()):
                        if value and key not in ["sample_xml"]:  # Skip XML and empty values
                            print(f"{key}: {value}")

            elif accession_type in ["assembly", "sequence", "protein"]:
                # Display assembly/sequence/protein metadata
                print(f"\n{accession_type.title()} Metadata Summary:")
                for key, value in sorted(metadata.items()):
                    if value and key not in ["ftp_path"]:  # Skip empty values and FTP path
                        print(f"{key}: {value}")

        else:
            logger.warning("No metadata retrieved")

    return 0


if __name__ == "__main__":
    sys.exit(main())
