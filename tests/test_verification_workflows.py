"""Keep ordinary CI numerical and expensive timing profiles explicitly manual."""

from pathlib import Path
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
    assert "--baseline" not in commands and "--profile release" not in commands
    assert data["permissions"] == {"contents": "read"}


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
