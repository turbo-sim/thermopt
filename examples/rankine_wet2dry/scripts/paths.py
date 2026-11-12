"""Shared locations; importing this module does not create files or run a case."""
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
STUDY = SCRIPTS.parent
PROJECT = STUDY.parents[1]
OUTPUT = STUDY / "output"
RESULTS = OUTPUT / "results"
LOGS = OUTPUT / "logs"
ARCHIVE = OUTPUT / "archive"
REFERENCE_INPUTS = OUTPUT / "reference_inputs"
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
        return OUTPUT / path
    return current
