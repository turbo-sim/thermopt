"""Verify real process-affinity control and inheritance outside numerical code.

Only a fresh child process is constrained. A grandchild observes the inherited
kernel mask, while the calling process retains its original processor set.
This is an explicit benchmark/CI check, with no JAX or thermodynamic imports.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if __package__ in (None, ""):
    sys.path.insert(0, str(ROOT))

from benchmark.run import cpu_affinity


def observe_worker(mask):
    """Constrain this disposable worker and observe its new child's mask."""
    before = cpu_affinity()
    constrained = cpu_affinity(mask)
    command = [
        sys.executable,
        "-c",
        "import json,sys; sys.path.insert(0,sys.argv[1]); "
        "from benchmark.run import cpu_affinity; "
        "print(json.dumps(cpu_affinity()))",
        str(ROOT),
    ]
    grandchild = subprocess.run(
        command, check=True, capture_output=True, text=True, timeout=30
    )
    inherited = json.loads(grandchild.stdout)
    assert int(constrained["actual_mask"], 16) == mask
    assert int(inherited["actual_mask"], 16) == mask
    return dict(before=before, constrained=constrained, grandchild=inherited)


def verify_inheritance(*, require_supported=False):
    """Check the actual OS through workers, preserving this process's mask."""
    before = cpu_affinity()
    if before["actual_mask"] is None:
        if require_supported:
            raise RuntimeError(
                "Real CPU-affinity observation is unavailable on this OS."
            )
        return dict(status="unsupported", parent_before=before)
    inherited = int(before["actual_mask"], 16)
    mask = inherited & -inherited
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker-mask",
        hex(mask),
    ]
    child = subprocess.run(
        command, check=True, capture_output=True, text=True, timeout=30
    )
    observed = json.loads(child.stdout)
    after = cpu_affinity()
    assert before == after, "The parent process's affinity changed."
    assert observed["before"]["actual_mask"] == before["actual_mask"]
    assert int(observed["constrained"]["actual_mask"], 16) == mask
    assert int(observed["grandchild"]["actual_mask"], 16) == mask
    return dict(
        status="passed",
        parent_before=before,
        requested_worker_mask=hex(mask),
        child=observed,
        parent_after=after,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-supported", action="store_true")
    parser.add_argument(
        "--worker-mask", type=lambda value: int(value, 0), help=argparse.SUPPRESS
    )
    args = parser.parse_args()
    if args.worker_mask is not None:
        print(json.dumps(observe_worker(args.worker_mask)))
        return
    if args.output is None:
        parser.error("--output is required")
    if args.output.exists():
        parser.error("Existing verification evidence is not overwritten")
    record = verify_inheritance(require_supported=args.require_supported)
    record.update(
        python=platform.python_version(),
        platform=platform.platform(),
        logical_cpus=os.cpu_count(),
        protocol_sha256=hashlib.sha256(
            (ROOT / "benchmark/run.py").read_bytes()
        ).hexdigest(),
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
