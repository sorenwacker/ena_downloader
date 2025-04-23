"""
Command-line interface for the ENA Downloader.
"""

import argparse
import logging
import os
import sys
from urllib.parse import urlparse

from ena_downloader.downloader import ENADownloader

logger = logging.getLogger(__name__)


def parse_args(args=None):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Download datasets from the European Nucleotide Archive (ENA)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("accessions", nargs="+", help="Study or project accession(s) to download (e.g., PRJNA123456)")

    parser.add_argument("-o", "--output-dir", default="data", help="Directory to save downloaded files")

    parser.add_argument(
        "--no-resume", action="store_true", help="Disable resume capability (download files from scratch)"
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    parser.add_argument("--list-only", action="store_true", help="Only list the files without downloading them")

    parser.add_argument("--clean-failed", action="store_true", help="Remove files that failed MD5 verification")

    parser.add_argument(
        "-m", "--max-files", type=int, default=None, help="Maximum number of files to download (randomly selected)"
    )

    parser.add_argument(
        "--max-per-accession",
        type=int,
        default=None,
        help="Maximum number of files to download per accession (randomly selected)",
    )

    return parser.parse_args(args)


def setup_logging(verbose: bool = False):
    """Set up logging configuration."""
    log_level = logging.DEBUG if verbose else logging.INFO

    # Configure root logger
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Suppress verbose logging from requests and urllib3
    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)


def display_files(dataset_info, accession):
    """Display available files for a given dataset."""
    logger.info("Files available for %s:", accession)

    for entry in dataset_info:
        if entry.get("fastq_ftp"):
            for url in entry["fastq_ftp"].split(";"):
                filename = os.path.basename(urlparse(url).path)
                logger.info("  %s (sequencing)", filename)

        if entry.get("submitted_ftp"):
            for url in entry["submitted_ftp"].split(";"):
                filename = os.path.basename(urlparse(url).path)
                logger.info("  %s (submitted)", filename)


def clean_failed_files(downloader):
    """Clean up files that failed MD5 verification."""
    downloaded_files = downloader._get_downloaded_files()
    failed_files = [
        path for path, info in downloaded_files.items() if info.get("status") in ("failed_md5", "failed_download")
    ]

    if failed_files:
        logger.info("Cleaning up %d failed files...", len(failed_files))
        for path in failed_files:
            if os.path.exists(path):
                try:
                    os.remove(path)
                    logger.info("Removed: %s", path)
                except Exception as e:
                    logger.error("Failed to remove %s: %s", path, e)
    else:
        logger.info("No failed files to clean up")


def main(argv=None):
    """Main function to run the ENA downloader from the command line."""
    args = parse_args(argv)

    # Set up logging
    setup_logging(args.verbose)

    # Initialize downloader with command-line arguments
    downloader = ENADownloader(output_dir=args.output_dir, resume=not args.no_resume, verbose=args.verbose)

    all_downloaded_files = []
    total_files_to_download = args.max_files

    # Process each accession
    for accession in args.accessions:
        try:
            logger.info("Processing dataset: %s", accession)

            # Get dataset information
            dataset_info = downloader.search_dataset(accession)

            if not dataset_info:
                logger.warning("No data found for accession: %s", accession)
                continue

            if args.list_only:
                # Just list the files without downloading
                display_files(dataset_info, accession)
            else:
                # Calculate max files for this accession
                max_files_this_accession = args.max_per_accession

                if total_files_to_download is not None:
                    # If global limit exists, make sure we don't exceed it
                    if max_files_this_accession is None:
                        max_files_this_accession = total_files_to_download
                    else:
                        max_files_this_accession = min(max_files_this_accession, total_files_to_download)

                # Download the dataset
                downloaded_files = downloader.download_dataset(accession, max_files=max_files_this_accession)
                all_downloaded_files.extend(downloaded_files)

                # Update the remaining files limit if applicable
                if total_files_to_download is not None:
                    total_files_to_download -= len(downloaded_files)
                    if total_files_to_download <= 0:
                        logger.info("Reached maximum number of files to download")
                        break

                if downloaded_files:
                    logger.info("Downloaded %d files for %s", len(downloaded_files), accession)
                else:
                    logger.warning("No files downloaded for %s", accession)

        except Exception as e:
            logger.error("Error processing dataset %s: %s", accession, e)
            if args.verbose:
                import traceback

                logger.debug(traceback.format_exc())

    # Print download statistics
    stats = downloader.get_download_stats()
    logger.info("Download Statistics:")
    logger.info("  Total files tracked: %d", stats["total"])
    logger.info("  Complete: %d", stats["complete"])
    logger.info("  Failed MD5: %d", stats["failed_md5"])
    logger.info("  Failed Download: %d", stats["failed_download"])
    logger.info("  Other status: %d", stats["other"])

    # Clean up failed files if requested
    if args.clean_failed:
        clean_failed_files(downloader)

    # Return success if we downloaded at least one file or just listed files
    if args.list_only or all_downloaded_files:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
