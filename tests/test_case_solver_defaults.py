"""Exercise the maintained YAML defaults through complete native solver runs."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
import thermopt as th
import yaml

from regression_support import EXAMPLE_CASES

ROOT = Path(__file__).resolve().parents[1]
MIGRATION = json.loads(
    (ROOT / "tests/data/case_solver_defaults_20261116.json").read_text()
)


@pytest.mark.parametrize("filename", sorted(MIGRATION["cases"]))
def test_default_solver_change_preserves_each_physical_formulation(filename):
    data = yaml.safe_load((ROOT / filename).read_text(encoding="utf-8"))
    digest = hashlib.sha256(
        json.dumps(data["problem_formulation"], sort_keys=True).encode()
    ).hexdigest()
    assert digest == MIGRATION["cases"][filename]["problem_formulation_sha256"]
    assert data["solver_options"]["library"] == "pysolver"
    assert data["solver_options"]["update_on"] == "iteration"


@pytest.mark.optimization
@pytest.mark.parametrize(
    "name", ["sco2_recuperated", "refrigeration_r134a", "orc_subcritical"]
)
def test_standard_case_defaults_complete_a_physical_controlled_slsqp_solve(
    name, tmp_path
):
    config = th.read_configuration_file(
        str(ROOT / "case_studies" / EXAMPLE_CASES[name])
    )
    config["solver_options"].update(
        callbacks={}, print_convergence=False, plot_convergence=False
    )
    source = tmp_path / "case.yaml"
    # Variable order affects finite-difference SLSQP numerics. Exercise the
    # launcher's actual order rather than safe_dump's alphabetical default.
    source.write_text(
        yaml.safe_dump(
            th.convert_numpy_to_python(config, precision=14), sort_keys=False
        )
    )
    manager = th.ThermodynamicCycleOptimization(str(source), out_dir=str(tmp_path))
    assert manager.problem.variable_names == list(
        config["problem_formulation"]["design_variables"]
    )
    manager.run_optimization()
    assert manager.solver.success, manager.solver.message
    result = manager.problem.fitness(manager.solver.x_final)
    assert np.all(np.isfinite(result))
    tolerance = float(config["solver_options"]["extra_options"]["feasibility_tol"])
    assert np.max(np.abs(manager.problem.c_eq), initial=0) <= tolerance
    assert np.max(manager.problem.c_ineq, initial=0) <= tolerance
    reference = json.loads(
        (ROOT / "tests/data/optimization_baselines.json").read_text()
    )["cases"][name]
    assert float(result[0]) == pytest.approx(
        reference["final"]["residuals"][0], rel=5e-5, abs=2e-6
    )
