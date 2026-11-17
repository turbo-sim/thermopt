"""Zero heat transfer must not derive a flow ratio from EOS roundoff."""

import numpy as np
import pytest
import jaxprop as jp

from thermopt.components import basic_components as basic
from thermopt.cycles import cycle_power_recuperated, cycle_refrigeration_recuperated
from regression_support import make_problem, read_baselines


@pytest.mark.parametrize(
    "name,point,module",
    [
        ("orc_transcritical", "initial", cycle_power_recuperated),
        ("heat_pump_co2_valve", "initial", cycle_refrigeration_recuperated),
        ("heat_pump_co2_valve", "nearby", cycle_refrigeration_recuperated),
        ("heat_pump_co2_turbine", "initial", cycle_refrigeration_recuperated),
        ("heat_pump_co2_turbine", "nearby", cycle_refrigeration_recuperated),
    ],
)
def test_zero_duty_recuperator_retains_its_unit_ratio_convention(
    name, point, module, tmp_path, monkeypatch
):
    records = []

    def observe(*args, **kwargs):
        result = basic.heat_exchanger(*args, **kwargs)
        records.append((result, args))
        return result

    monkeypatch.setattr(module, "heat_exchanger", observe)
    expected = read_baselines()["cases"][name]
    design = expected["x0"] if point == "initial" else expected["nearby_x"]
    problem = make_problem(name, tmp_path)
    problem.fitness(design)
    recuperator = problem.cycle_data["components"]["recuperator"]
    assert any(result is recuperator for result, _ in records)
    # Earlier EOS flashes can perturb even the supplied endpoint enthalpies.
    # The bypass is defined by the cycle, and its independently assigned stream
    # flows establish the ratio without dividing tiny enthalpy differences.
    if name == "orc_transcritical":
        assert problem.x0_dict["recuperator_effectiveness"] == 0.0
    else:
        assert problem.cycle_topology == "refrigeration_simple"
    hot_flow = recuperator["hot_side"]["mass_flow"]
    cold_flow = recuperator["cold_side"]["mass_flow"]
    assert hot_flow == cold_flow and hot_flow > 0.0
    assert recuperator["mass_flow_ratio"] == 1.0


def test_small_nonzero_recuperation_remains_active(tmp_path):
    problem = make_problem("sco2_recuperated", tmp_path)
    design = problem.x0.copy()
    design[problem.variable_names.index("recuperator_effectiveness")] = 1e-8
    problem.fitness(design)
    recuperator = problem.cycle_data["components"]["recuperator"]
    assert recuperator["q_hot_side"] > 0.0
    assert recuperator["q_cold_side"] > 0.0
    assert np.isfinite(recuperator["mass_flow_ratio"])
    assert recuperator["mass_flow_ratio"] == pytest.approx(
        recuperator["q_cold_side"] / recuperator["q_hot_side"], rel=1e-12
    )


def test_generic_exchanger_preserves_a_small_nonzero_nonunit_duty_ratio():
    fluid = jp.Fluid("Air", method="coolprop", backend="HEOS")
    hot_in, cold_in = 500e3, 300e3
    hot_out, cold_out = hot_in - 1e-2, cold_in + 2e-2
    exchanger = basic.heat_exchanger(
        fluid,
        hot_in,
        hot_out,
        2e5,
        2e5,
        fluid,
        cold_in,
        cold_out,
        1e5,
        1e5,
        num_steps=3,
    )
    assert exchanger["q_hot_side"] > 0.0
    assert exchanger["q_cold_side"] > 0.0
    assert exchanger["mass_flow_ratio"] == pytest.approx(2.0, rel=2e-7)
