"""Keep ordinary CI numerical and expensive timing profiles explicitly manual."""

from pathlib import Path
import json
import os
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


def workflow(name):
    # GitHub uses YAML 1.2: preserve its `on` key instead of YAML 1.1 booleans.
    return yaml.load(
        (ROOT / ".github/workflows" / name).read_text(), Loader=yaml.BaseLoader
    )


def test_correctness_matrix_covers_supported_interpreters_and_both_platforms():
    data = workflow("tests.yaml")
    matrix = data["jobs"]["pytest"]["strategy"]["matrix"]
    assert set(matrix["python-version"]) == {"3.11", "3.13"}
    assert set(matrix["os"]) == {"ubuntu-latest", "windows-latest"}
    steps = data["jobs"]["pytest"]["steps"]
    commands = "\n".join(step.get("run", "") for step in steps)
    assert "poetry check --lock" in commands and "--only main" in commands
    assert "--junitxml" in commands and "--profile smoke" in commands
    assert "verify_affinity.py --require-supported" in commands
    assert "--baseline" not in commands and "--profile release" not in commands
    assert data["permissions"] == {"contents": "read"}


def test_pypi_publication_requires_a_version_tag_and_a_supported_lock_tool():
    data = workflow("publish_to_pypi.yaml")
    assert set(data["on"]) == {"push"}
    assert data["on"]["push"] == {"tags": ["v*.*.*"]}
    commands = "\n".join(
        step.get("run", "") for step in data["jobs"]["publish_to_pypi"]["steps"]
    )
    assert "poetry==2.2.1" in commands and "poetry check --lock" in commands
    assert data["permissions"] == {"contents": "read"}


@pytest.mark.skipif(
    os.name != "nt" and not hasattr(os, "sched_getaffinity"),
    reason="This OS does not provide the supported affinity API",
)
def test_real_affinity_is_inherited_by_workers_without_changing_the_parent(tmp_path):
    output = tmp_path / "affinity.json"
    process = subprocess.run(
        [
            sys.executable,
            str(ROOT / "benchmark/verify_affinity.py"),
            "--require-supported",
            "--output",
            str(output),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=45,
    )
    assert process.returncode == 0, process.stdout + process.stderr
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["status"] == "passed"
    assert result["parent_before"] == result["parent_after"]
    requested = result["requested_worker_mask"]
    assert result["child"]["constrained"]["actual_mask"] == requested
    assert result["child"]["grandchild"]["actual_mask"] == requested


def test_benchmark_workflow_is_opt_in_and_retains_measurements():
    data = workflow("benchmarks.yaml")
    assert set(data["on"]) == {"workflow_dispatch"}
    profile = data["on"]["workflow_dispatch"]["inputs"]["profile"]
    assert set(profile["options"]) == {"smoke", "phase", "release"}
    steps = data["jobs"]["benchmark"]["steps"]
    artifacts = [
        step
        for step in steps
        if step.get("uses", "").startswith("actions/upload-artifact@")
    ]
    assert artifacts and all(step["if"] == "always()" for step in artifacts)
    assert data["permissions"] == {"contents": "read"}
