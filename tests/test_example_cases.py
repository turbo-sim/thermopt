"""Cycle regressions replace running interactive example scripts in-place."""

import copy

import numpy as np
import pytest

from thermopt.optimize_cycle import CYCLE_TOPOLOGIES

from regression_support import (EXAMPLE_CASES, assert_snapshot, make_problem,
                                read_baselines, snapshot)


@pytest.mark.parametrize("name", list(EXAMPLE_CASES))
@pytest.mark.parametrize("point", ["initial", "nearby"])
def test_cycle_numerical_regression(name, point, tmp_path):
    """Freeze residual ordering, states, energy metrics and complete HX curves."""
    expected = read_baselines()["cases"][name]
    assert "recording_error" not in expected, f"Missing reference for {name}"
    problem = make_problem(name, tmp_path)
    assert problem.variable_names == expected["variable_names"]
    assert problem.cycle_topology == expected["topology"]
    np.testing.assert_allclose(problem.x0, expected["x0"], rtol=1e-12, atol=1e-8)
    x = expected["x0"] if point == "initial" else expected["nearby_x"]
    actual = snapshot(problem, x)
    assert_snapshot(actual, expected[point])
    values = np.asarray(actual["residuals"])
    assert values.shape == (1 + problem.get_nec() + problem.get_nic(),)
    np.testing.assert_array_equal(values[1:1 + problem.get_nec()], problem.c_eq)
    np.testing.assert_array_equal(values[1 + problem.get_nec():], problem.c_ineq)
    assert len(problem.constraint_data) == values.size - 1


@pytest.mark.parametrize("name", list(EXAMPLE_CASES))
def test_cycle_component_and_system_energy_balances(name, tmp_path):
    """Physical checks complement characterization of the existing outputs."""
    problem = make_problem(name, tmp_path)
    problem.fitness(problem.x0)
    energy = problem.cycle_data["energy_analysis"]
    scale = max(1.0, *(abs(float(v)) for k, v in energy.items() if "power" in k or "heat_flow" in k))
    for key, value in energy.items():
        if "energy_balance" in key:
            # This report mixes source/sink pump power into the working-loop
            # residual for recompression only; the separate xfail tracks it.
            if name == "sco2_recompression":
                value -= energy["heat_source_pump_power"] + energy["heat_sink_pump_power"]
            assert abs(value) <= 1e-7 * scale, (name, key, value)
    for component_name, component in problem.cycle_data["components"].items():
        if component["type"] == "heat_exchanger":
            hot, cold = component["hot_side"], component["cold_side"]
            hot_duty = hot["mass_flow"] * (hot["state_in"].h - hot["state_out"].h)
            cold_duty = cold["mass_flow"] * (cold["state_out"].h - cold["state_in"].h)
            np.testing.assert_allclose(hot_duty, cold_duty, rtol=2e-8, atol=0.05,
                                       err_msg=f"{name}/{component_name}")
            np.testing.assert_allclose(component["temperature_difference"],
                                       hot["states"].T - cold["states"].T, atol=1e-10)
        elif component["type"] in ["compressor", "expander"]:
            sign = 1 if component["type"] == "compressor" else -1
            work = sign * (component["state_out"].h - component["state_in"].h)
            np.testing.assert_allclose(component["specific_work"], work, rtol=1e-10, atol=1e-7)
            np.testing.assert_allclose(component["power"], work * component["mass_flow"], rtol=1e-10)


@pytest.mark.parametrize("name", ["sco2_recuperated", "refrigeration_r134a", "ptes", "ptes_turbo"])
def test_repeated_evaluation_has_no_cross_design_state_leak(name, tmp_path):
    expected = read_baselines()["cases"][name]
    problem = make_problem(name, tmp_path)
    initial = snapshot(problem, expected["x0"])
    snapshot(problem, expected["nearby_x"])
    repeated = snapshot(problem, expected["x0"])
    assert_snapshot(repeated, initial)
    np.testing.assert_array_equal(list(problem.x0_dict.values()), expected["x0"])


@pytest.mark.parametrize("canonical,alias,case", [
    ("simple", "power_simple", "orc_subcritical"),
    ("recuperated", "power_recuperated", "sco2_recuperated"),
    ("recompression", "power_split_compression", "sco2_recompression"),
    ("recompression", "split_compression", "sco2_recompression"),
])
def test_public_topology_aliases_preserve_outputs(canonical, alias, case, tmp_path):
    assert CYCLE_TOPOLOGIES[canonical] is CYCLE_TOPOLOGIES[alias]
    problem = make_problem(case, tmp_path)
    original = snapshot(problem, problem.x0)
    problem.cycle_topology = alias
    assert_snapshot(snapshot(problem, problem.x0), original)


def test_cycle_evaluator_preserves_caller_inputs(tmp_path):
    problem = make_problem("sco2_recuperated", tmp_path)
    variables = dict(zip(problem.variable_names, problem.x0))
    parameters = copy.deepcopy(problem.fixed_parameters)
    variables_before, parameters_before = copy.deepcopy((variables, parameters))
    CYCLE_TOPOLOGIES[problem.cycle_topology](variables, parameters, problem.constraints, problem.objective_function)
    np.testing.assert_equal(variables, variables_before)
    np.testing.assert_equal(parameters, parameters_before)


def test_unknown_topology_is_reported(tmp_path):
    problem = make_problem("sco2_recuperated", tmp_path)
    problem.cycle_topology = "does-not-exist"
    with pytest.raises(ValueError, match="Invalid cycle topology"):
        problem.fitness(problem.x0)


def test_wet_expansion_fixture_exercises_two_phase_inlet_and_outlet(tmp_path):
    problem = make_problem("orc_wet_expansion", tmp_path)
    problem.fitness(problem.x0)
    expander = problem.cycle_data["components"]["expander"]
    assert expander["state_in"].is_two_phase == 1
    assert expander["state_out"].is_two_phase == 1
    assert expander["state_in"].Q == pytest.approx(0.7, abs=1e-10)
    assert 0 < expander["state_in"].Q < expander["state_out"].Q < 1


@pytest.mark.xfail(strict=True, raises=AssertionError, reason="Recompression energy_balance includes auxiliary pump work outside the working-fluid heat duties.")
def test_recompression_reported_energy_balance_closes(tmp_path):
    problem = make_problem("sco2_recompression", tmp_path)
    problem.fitness(problem.x0)
    assert abs(problem.cycle_data["energy_analysis"]["energy_balance"]) < 0.1
