"""Run the suite from any directory, forwarding arguments to pytest."""

from pathlib import Path
import sys

import pytest


if __name__ == "__main__":
    raise SystemExit(pytest.main([str(Path(__file__).resolve().parent), *sys.argv[1:]]))
