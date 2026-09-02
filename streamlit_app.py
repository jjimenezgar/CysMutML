"""Entry point for the CysMutML Streamlit application."""

from pathlib import Path
import sys

# Streamlit Cloud runs this root-level file directly. The package uses the
# standard src layout, so make that package root explicit for Cloud and local launches.
PACKAGE_ROOT = Path(__file__).resolve().parent / "src"
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from cysmutml.webapp import main  # noqa: E402

main()
