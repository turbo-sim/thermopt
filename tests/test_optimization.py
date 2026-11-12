"""Complete SLSQP solves, checked by physics/feasibility rather than iterations."""

import json

import numpy as np
import pytest

from regression_support import DATA, OPTIMIZATION_CASES, make_optimization


@pytest.mark.optimization
@pytest.mark.slow
@pytest.mark.parametrize("name", OPTIMIZATION_CASES)
def test_complete_cycle_optimization_matches_reference(name, tmp_path):
    reference = json.loads((DATA / "optimization_baselines.json").read_text(encoding="utf-8"))["cases"][name]
    assert reference["success"], "Reference must be a successful complete solve"
    manager = make_optimization(name, tmp_path)
    problem = manager.problem
    initial = problem.x0.copy()
    initial_residuals = problem.fitness(initial).copy()
    manager.run_optimization()
    solver = manager.solver
    assert solver.success, solver.message
    assert 0 < solver.iteration_count < 200
    assert solver.func_count_tot > 1  # A real solve, not a mocked or saved result.

    point = np.asarray(solver.x_final, dtype=float)
    assert point.shape == initial.shape
    assert np.all(np.isfinite(point))
    lower, upper = map(np.asarray, problem.get_bounds())
    normalized_violation = np.maximum(lower - point, point - upper) / (upper - lower)
    assert np.max(normalized_violation) < 1e-7

    # Reevaluate the solver's returned physical design; do not trust cached
    # objective/feasibility from a previous finite-difference probe.
    final = np.asarray(problem.fitness(point))
    assert np.all(np.isfinite(final))
    assert final.shape == initial_residuals.shape
    assert np.max(np.abs(problem.c_eq), initial=0.0) <= 2e-6
    assert np.max(problem.c_ineq, initial=0.0) <= 2e-6
    expected = reference["final"]
    assert final[0] == pytest.approx(expected["residuals"][0], rel=5e-5, abs=2e-6)
    # These are objective-related physical metrics; inactive design variables
    # can have multiple optima, so exact x/iteration-history equality is avoided.
    metric = "COP_refrigeration" if name == "refrigeration_r134a" else "cycle_efficiency"
    actual_metric = float(problem.cycle_data["energy_analysis"][metric])
    assert actual_metric == pytest.approx(expected["energy_analysis"][metric], rel=5e-5, abs=2e-6)
    assert actual_metric > 0
    assert not np.allclose(point, initial, rtol=1e-9, atol=1e-9)

    # Operating states are protected by tight independent regressions at the
    # saved optimum, without assuming every optimizer picks the same free point.
    saved = np.asarray(problem.fitness(reference["x_final"]))
    np.testing.assert_allclose(saved, expected["residuals"], rtol=2e-6, atol=2e-6)
    for component_name in ["compressor", "expander"]:
        component = problem.cycle_data["components"][component_name]
        for endpoint in ["state_in", "state_out"]:
            for prop, atol in [("p", 0.1), ("h", 0.02), ("T", 2e-5), ("s", 2e-4), ("rho", 1e-6)]:
                assert float(component[endpoint][prop]) == pytest.approx(
                    expected["components"][component_name][endpoint][prop], rel=2e-7, abs=atol)
