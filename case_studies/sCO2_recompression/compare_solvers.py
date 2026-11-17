"""Case-local launcher; comparison implementation lives under benchmark/."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchmark.recompression_solvers import main

if __name__ == "__main__":
    main()
