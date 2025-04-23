"""
Utility functions for the ENA Downloader.
"""

import hashlib
import logging
import os

logger = logging.getLogger(__name__)


def calculate_md5(file_path: str, chunk_size: int = 4096) -> str:
    """
    Calculate MD5 checksum of a file.

    Args:
        file_path: Path to the file
        chunk_size: Size of chunks to read

    Returns:
        MD5 checksum as a hex string
    """
    md5_hash = hashlib.md5()

    try:
        with open(file_path, "rb") as f:
            # Read the file in chunks
            for chunk in iter(lambda: f.read(chunk_size), b""):
                md5_hash.update(chunk)

        return md5_hash.hexdigest()
    except Exception as e:
        logger.error("Error calculating MD5 for %s: %s", file_path, e)
        return ""


def get_file_size(file_path: str) -> int:
    """
    Get the size of a file in bytes.

    Args:
        file_path: Path to the file

    Returns:
        Size of the file in bytes
    """
    try:
        return os.path.getsize(file_path)
    except Exception as e:
        logger.error("Error getting file size for %s: %s", file_path, e)
        return 0


def human_readable_size(size_bytes: int) -> str:
    """
    Convert size in bytes to a human-readable format.

    Args:
        size_bytes: Size in bytes

    Returns:
        Human-readable string representation of the size
    """
    if size_bytes == 0:
        return "0B"

    size_names = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = 0
    while size_bytes >= 1024 and i < len(size_names) - 1:
        size_bytes /= 1024.0
        i += 1

    return f"{size_bytes:.2f} {size_names[i]}"


def scan_directory(directory: str) -> list[dict[str, str]]:
    """
    Scan a directory for files and calculate their MD5 checksums.

    Args:
        directory: Directory to scan

    Returns:
        List of dictionaries with file information
    """
    file_info = []

    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            file_size = get_file_size(file_path)

            info = {
                "path": file_path,
                "size": file_size,
                "size_human": human_readable_size(file_size),
                "md5": calculate_md5(file_path),
            }

            file_info.append(info)

    return file_info


def parse_accession(accession: str) -> tuple[str, str]:
    """
    Parse an accession string to determine its type.

    Args:
        accession: Accession string (e.g., PRJEB12345, SRR123456)

    Returns:
        Tuple of (accession_type, accession)
    """
    prefix = accession[:3]

    if prefix == "PRJ":
        return "project", accession
    if prefix == "SRP":
        return "study", accession
    if prefix == "SRR":
        return "run", accession
    if prefix == "SRS":
        return "sample", accession
    if prefix == "ERR":
        return "run", accession
    if prefix == "ERP":
        return "study", accession
    if prefix == "DRR":
        return "run", accession
    if prefix == "DRP":
        return "study", accession
    return "unknown", accession
