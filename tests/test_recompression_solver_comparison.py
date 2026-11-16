"""Comparison setup must not rewrite physical inputs or leak solver controls."""

import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import thermopt as th

from regression_support import load_case

RUNNER = (
    Path(__file__).resolve().parents[1]
    / "case_studies/sCO2_recompression/compare_solvers.py"
)
spec = importlib.util.spec_from_file_location("recompression_comparison", RUNNER)
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_kkt_export_accepts_native_integer_keys_and_preserves_gaps(tmp_path):
    path = tmp_path / "kkt.json"
    comparison.save_json(
        path, {"multipliers": {np.int64(2): np.float64(3)}, "unavailable": np.nan}
    )
    assert json.loads(path.read_text()) == {
        "multipliers": {"2": 3.0},
        "unavailable": None,
    }


@pytest.mark.parametrize("backend", ["coolprop", "bicubic"])
def test_grid_override_preserves_yaml_recipe_and_physical_inputs(tmp_path, backend):
    source = RUNNER.parent / "case_sCO2_recompression.yaml"
    before = source.read_bytes()
    configured = th.read_configuration_file(str(source))
    original = copy.deepcopy(configured)
    actual = comparison.comparison_configuration(
        configured, tmp_path / "tables", backend
    )
    assert source.read_bytes() == before
    np.testing.assert_equal(configured, original)
    fluid = actual["problem_formulation"]["fixed_parameters"]["working_fluid"]
    assert fluid["method"] == backend
    if backend == "bicubic":
        assert fluid["table_options"]["N_h"] == fluid["table_options"]["N_p"] == 128
    else:
        assert "table_options" not in fluid
    assert (
        original["problem_formulation"]["fixed_parameters"]["working_fluid"][
            "table_options"
        ]["N_h"]
        == 200
    )
    for key in ("design_variables", "constraints", "objective_function"):
        np.testing.assert_equal(
            actual["problem_formulation"][key], original["problem_formulation"][key]
        )
    for role in ("heating_fluid", "cooling_fluid"):
        np.testing.assert_equal(
            actual["problem_formulation"]["fixed_parameters"][role],
            original["problem_formulation"]["fixed_parameters"][role],
        )


@pytest.mark.parametrize(
    "library,method", [("pysolver", "slsqp"), ("pounce", "ipm"), ("pounce", "sqp")]
)
def test_each_solver_gets_common_numerics_and_supported_controls(library, method):
    config = th.read_configuration_file(
        str(RUNNER.parent / "case_sCO2_recompression.yaml")
    )
    actual = comparison.solver_options(config, library, method, None, False)
    assert type(actual["max_iterations"]) is int
    for key in (
        "problem_scale",
        "derivative_method",
        "derivative_abs_step",
        "tolerance",
    ):
        assert actual[key] == config["solver_options"][key]
    assert "callbacks" not in actual
    if library == "pounce":
        assert (
            not {"first_step_fraction", "max_step_fraction", "step_growth_factor"}
            & actual["extra_options"].keys()
        )
    else:
        assert actual["extra_options"]["first_step_fraction"] == 0.2


def test_manager_accepts_yaml_integer_limit_with_controlled_slsqp(tmp_path):
    config = load_case("sco2_recuperated")
    config["solver_options"].update(
        library="pysolver",
        max_iterations=np.int64(1),
        callbacks={},
        print_convergence=False,
        plot_convergence=False,
        extra_options=dict(first_step_fraction=0.01, optimality_tol=1e-4),
    )
    source = tmp_path / "controlled.yaml"
    import yaml

    source.write_text(
        yaml.safe_dump(th.convert_numpy_to_python(config, precision=14)),
        encoding="utf-8",
    )
    manager = th.ThermodynamicCycleOptimization(str(source), out_dir=str(tmp_path))
    manager.solver.solve(manager.problem.x0)
    assert manager.solver.result is not None
    assert manager.solver.iteration_count <= 1
