"""Historical entry point; maintained native workloads live under benchmark/.

Dated reports retain their original source hashes. New runs use benchmark
protocol 3 and never normalize obsolete JaxProp constructor settings.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark.native_backends import main

if __name__ == "__main__":
    main()
