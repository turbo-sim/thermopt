"""Headless, isolated test execution; do not write results into case_studies/."""

import os
from pathlib import Path
import shutil
import sys
from uuid import uuid4

os.environ["MPLBACKEND"] = "Agg"
os.environ["DISABLE_PLOTS"] = "1"

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import pytest


@pytest.fixture(scope="session")
def test_workspace():
    """Ordinary directory permissions also work in restricted Windows sessions.

    Python/pytest's private mode-0700 temp directories can be inaccessible to
    a managed Windows process. Keep artifacts isolated without relying on them.
    """
    parent = ROOT / ".test_artifacts"
    parent.mkdir(exist_ok=True)
    directory = parent / uuid4().hex
    directory.mkdir()
    yield directory
    # Only delete this run's explicitly verified child, never an arbitrary path.
    resolved = directory.resolve()
    if resolved.parent != parent.resolve():
        raise RuntimeError(f"Unexpected test artifact path: {resolved}")
    shutil.rmtree(resolved)


@pytest.fixture
def tmp_path(test_workspace):
    """Provide pytest's familiar isolated Path without mode-0700 directories."""
    directory = test_workspace / uuid4().hex
    directory.mkdir()
    return directory


@pytest.fixture(autouse=True)
def close_figures():
    """Close figures even if an assertion fails; restore global style changes."""
    with matplotlib.rc_context():
        yield
    plt.close("all")
