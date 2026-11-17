"""Shared locations; importing this module does not create files or run a case."""
from pathlib import Path
import os

SCRIPTS = Path(__file__).resolve().parent
STUDY = SCRIPTS.parent
PROJECT = STUDY.parents[1]
LEGACY_OUTPUT = STUDY / "output"
# New generated artifacts stay together. Existing research records remain at
# their original locations and are never moved or relabeled as current solves.
OUTPUT = Path(os.getenv("THERMOPT_RESULTS_DIR", STUDY / "results"))
RESULTS = OUTPUT / "studies"
LOGS = OUTPUT / "logs"
ARCHIVE = LEGACY_OUTPUT / "archive"
REFERENCE_INPUTS = LEGACY_OUTPUT / "reference_inputs"
WEB_ASSETS = SCRIPTS / "web_assets"
HTML_PAGES = OUTPUT / "report"


def resolve_record_path(value):
    """Resolve paths recorded before or after the directory reorganization.

    Original native solver records remain byte-for-byte unchanged. Their recorded
    initializer paths used results/, archive/ or reference_inputs/ at the root.
    """
    path = Path(value)
    if path.is_absolute():
        return path
    current = STUDY / path
    if current.exists():
        return current
    if path.parts and path.parts[0] in ("results", "archive", "reference_inputs", "report"):
        historical = LEGACY_OUTPUT / path
        return historical if historical.exists() else OUTPUT / path
    return current
