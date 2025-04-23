"""
Pytest configuration file for ENA Downloader tests.
"""

import sys
from pathlib import Path

# Add the src directory to the Python path
src_dir = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_dir))
