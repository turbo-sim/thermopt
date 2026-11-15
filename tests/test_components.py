"""Component contracts, physical balances and small numerical regressions.

Reference values were captured with jaxprop 0.5.7 / HEOS (CoolProp 7.2.0).
Only the public jaxprop property interface is used here. Tolerances permit
ordinary property-library roundoff while detecting changes in model behavior.
"""

import numpy as np
import pytest
import jaxprop as props

from thermopt.components.basic_components import (
    compression_process,
    compute_component_energy_flows,
    expansion_process,
    heat_exchanger,
    heat_transfer_process,
    isenthalpic_valve,
)
from thermopt.components.turbomachinery_nondimensional import (
    compressor_defaults,
    turbine_defaults,
    friction_factor,
)


@pytest.fixture(scope="module")
def air():
    return props.Fluid("Air", backend="HEOS")


def assert_finite_state(state):
    for name in ("p", "T", "h", "s", "rho"):
        assert np.all(np.isfinite(np.asarray(state[name]))), name
    assert np.all(np.asarray(state.T) > 0)
    assert np.all(np.asarray(state.p) > 0)
    assert np.all(np.asarray(state.rho) > 0)


def make_exchanger(air, counter_current=True, zero_duty=False):
    hot_in = air.get_state(props.PT_INPUTS, 2e5, 650.0)
    hot_out = air.get_state(props.PT_INPUTS, 1.9e5, 400.0)
    cold_in = air.get_state(props.PT_INPUTS, 1e5, 300.0)
    cold_out = air.get_state(props.PT_INPUTS, 0.95e5, 350.0)
    return heat_exchanger(
        air, hot_in.h, hot_in.h if zero_duty else hot_out.h, hot_in.p, hot_in.p if zero_duty else hot_out.p,
        air, cold_in.h, cold_in.h if zero_duty else cold_out.h, cold_in.p, cold_in.p if zero_duty else cold_out.p,
        num_steps=9, counter_current=counter_current,
    )


def test_heat_transfer_process_preserves_linear_pressure_enthalpy_path(air):
    inlet = air.get_state(props.PT_INPUTS, 2e5, 300.0)
    outlet = air.get_state(props.PT_INPUTS, 1.9e5, 500.0)
    path = heat_transfer_process(air, inlet.h, inlet.p, outlet.h, outlet.p, num_steps=7)

    assert path["fluid_name"] == "Air"
    assert_finite_state(path["states"])
    np.testing.assert_allclose(path["states"].h, np.linspace(inlet.h, outlet.h, 7))
    np.testing.assert_allclose(path["states"].p, np.linspace(inlet.p, outlet.p, 7))
    np.testing.assert_allclose([path["state_in"].T, path["state_out"].T], [300, 500], atol=1e-7)
    assert np.all(np.diff(path["states"].T) > 0)
    assert np.isnan(path["mass_flow"])


@pytest.mark.parametrize("counter_current", [False, True], ids=["cocurrent", "countercurrent"])
def test_exchanger_temperature_mapping_and_energy_ratio(air, counter_current):
    hx = make_exchanger(air, counter_current)
    hot, cold = hx["hot_side"], hx["cold_side"]
    assert_finite_state(hot["states"])
    assert_finite_state(cold["states"])
    # Physical inlet/outlet states retain the stream direction; only the hot
    # plotting/pinch array reverses for counter-current heat exchange.
    np.testing.assert_allclose([hot["state_in"].T, hot["state_out"].T], [650, 400], atol=1e-7)
    np.testing.assert_allclose(np.asarray(hot["states"].T)[[0, -1]], [400, 650] if counter_current else [650, 400], atol=1e-7)
    np.testing.assert_allclose(np.asarray(cold["states"].T)[[0, -1]], [300, 350], atol=1e-7)
    expected_gaps = [100, 300] if counter_current else [350, 50]
    np.testing.assert_allclose(np.asarray(hx["temperature_difference"])[[0, -1]], expected_gaps, atol=1e-7)
    assert np.all(hx["temperature_difference"] > 0)
    assert hx["q_hot_side"] > hx["q_cold_side"] > 0
    assert hx["mass_flow_ratio"] * hx["q_hot_side"] == pytest.approx(hx["q_cold_side"])


def test_zero_duty_exchanger_is_finite_and_has_unit_mass_ratio(air):
    hx = make_exchanger(air, zero_duty=True)
    assert hx["mass_flow_ratio"] == 1.0
    assert hx["q_hot_side"] == 0.0
    assert hx["q_cold_side"] == 0.0
    assert np.all(np.isfinite(hx["temperature_difference"]))


@pytest.mark.parametrize(
    "process,p_in,t_in,p_out,temperature,work",
    [
        (compression_process, 1e5, 300, 2e5, 377.12061335751343, 77648.19549879926),
        (expansion_process, 2e5, 600, 1e5, 511.23376852158447, 92425.66103170137),
    ],
    ids=["compressor", "expander"],
)
def test_isentropic_efficiency_and_air_outlet_regression(air, process, p_in, t_in, p_out, temperature, work):
    inlet = air.get_state(props.PT_INPUTS, p_in, t_in)
    component = process(air, inlet.h, inlet.p, p_out, efficiency=0.85)
    outlet = component["state_out"]
    assert_finite_state(outlet)
    assert outlet.p == pytest.approx(p_out)
    assert outlet.s > inlet.s  # Irreversible adiabatic process.
    assert component["specific_work"] > 0
    assert outlet.T == pytest.approx(temperature, rel=2e-6)
    assert component["specific_work"] == pytest.approx(work, rel=2e-6)
    if process is compression_process:
        measured_efficiency = component["isentropic_work"] / component["specific_work"]
    else:
        measured_efficiency = component["specific_work"] / component["isentropic_work"]
    assert measured_efficiency == pytest.approx(0.85, abs=1e-10)


@pytest.mark.parametrize(
    "process,p_in,t_in,p_out,temperature,work",
    [
        (compression_process, 9e6, 320, 18e6, 368.536578902871, 27293.022517137753),
        (expansion_process, 18e6, 650, 9e6, 580.0251393282734, 66631.70203815226),
    ],
    ids=["dense_CO2_compressor", "CO2_expander"],
)
def test_real_fluid_machinery_regression(process, p_in, t_in, p_out, temperature, work):
    fluid = props.Fluid("CO2", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, p_in, t_in)
    component = process(fluid, inlet.h, inlet.p, p_out, efficiency=0.85)
    assert_finite_state(component["state_out"])
    assert component["state_out"].T == pytest.approx(temperature, rel=2e-6)
    assert component["specific_work"] == pytest.approx(work, rel=2e-6)


@pytest.mark.parametrize(
    "process,p_in,t_in,p_out,temperature,work",
    [
        (compression_process, 1e5, 300, 2e5, 378.519037347315, 79064.23908188386),
        (expansion_process, 2e5, 600, 1e5, 509.9713671845537, 93728.35196805175),
    ],
    ids=["compressor", "expander"],
)
def test_polytropic_path_and_specific_work_regression(air, process, p_in, t_in, p_out, temperature, work):
    inlet = air.get_state(props.PT_INPUTS, p_in, t_in)
    component = process(air, inlet.h, inlet.p, p_out, efficiency=0.85,
                        efficiency_type="polytropic", num_steps=9)
    path = component["states"]
    assert_finite_state(path)
    assert np.asarray(path.p).shape == (9,)
    np.testing.assert_allclose(path.p, np.linspace(p_in, p_out, 9))
    assert np.all(np.diff(path.s) > 0)
    assert component["state_out"].T == pytest.approx(temperature, rel=2e-6)
    assert component["specific_work"] == pytest.approx(work, rel=2e-6)
    assert component["data_out"]["polytropic_efficiency"] == 0.85


@pytest.mark.parametrize("process,p_in,t_in,p_out", [
    (compression_process, 1e5, 300, 2e5),
    (expansion_process, 2e5, 600, 1e5),
], ids=["compressor", "expander"])
def test_reversible_limit_of_both_efficiency_models(air, process, p_in, t_in, p_out):
    inlet = air.get_state(props.PT_INPUTS, p_in, t_in)
    isentropic = process(air, inlet.h, inlet.p, p_out, efficiency=1.0)
    polytropic = process(air, inlet.h, inlet.p, p_out, efficiency=1.0,
                         efficiency_type="polytropic", num_steps=5)
    assert isentropic["state_out"].s == pytest.approx(inlet.s, abs=1e-7)
    assert polytropic["state_out"].T == pytest.approx(isentropic["state_out"].T, rel=1e-4)
    assert polytropic["specific_work"] == pytest.approx(isentropic["specific_work"], rel=1e-3)


@pytest.mark.parametrize("process,p_in,t_in,p_out", [
    (compression_process, 1e5, 300, 2e5),
    (expansion_process, 2e5, 600, 1e5),
], ids=["compressor", "expander"])
def test_unknown_efficiency_model_raises(air, process, p_in, t_in, p_out):
    inlet = air.get_state(props.PT_INPUTS, p_in, t_in)
    with pytest.raises(ValueError, match="efficiency_type"):
        process(air, inlet.h, p_in, p_out, 0.85, efficiency_type="unsupported")


@pytest.mark.parametrize("name,p_in,t_in,p_out", [
    ("Air", 2e5, 300, 1e5),
    ("R134a", 1e6, 300, 2e5),
], ids=["gas", "flashing_refrigerant"])
def test_valve_is_isenthalpic_and_irreversible(name, p_in, t_in, p_out):
    fluid = props.Fluid(name, backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, p_in, t_in)
    path = isenthalpic_valve(inlet, p_out, fluid, N=7)
    assert_finite_state(path)
    np.testing.assert_allclose(path.h, np.full(7, inlet.h), rtol=1e-10, atol=1e-5)
    np.testing.assert_allclose(path.p, np.linspace(p_in, p_out, 7))
    assert np.all(np.diff(path.s) > 0)
    if name == "R134a":
        assert path.is_two_phase[-1] == 1
        assert 0 < path.Q[-1] < 1


def test_component_energy_flows_conserve_heat_and_scale_shaft_work(air):
    hx = make_exchanger(air)
    hx["hot_side"]["mass_flow"] = 2.0
    hx["cold_side"]["mass_flow"] = 2.0 / hx["mass_flow_ratio"]
    inlet = air.get_state(props.PT_INPUTS, 1e5, 300.0)
    compressor = compression_process(air, inlet.h, inlet.p, 2e5, 0.85)
    compressor["mass_flow"] = 3.0
    inlet_turbine = air.get_state(props.PT_INPUTS, 2e5, 600.0)
    turbine = expansion_process(air, inlet_turbine.h, inlet_turbine.p, 1e5, 0.85)
    turbine["mass_flow"] = 3.0
    assert compute_component_energy_flows({"hx": hx, "compressor": compressor, "turbine": turbine}) is None
    assert hx["heat_flow"] == pytest.approx(2.0 * hx["q_hot_side"])
    assert hx["heat_flow_bis"] == pytest.approx(hx["heat_flow"])
    assert hx["heat_balance"] == pytest.approx(0.0, abs=1e-7)
    assert hx["power"] == 0.0
    for machine in (compressor, turbine):
        assert machine["power"] == pytest.approx(3.0 * machine["specific_work"])
        assert machine["heat_flow_rate"] == 0.0


@pytest.mark.xfail(strict=True, raises=AssertionError, reason="Existing co-current energy accounting uses hot heat_flow[0], which is zero.")
def test_cocurrent_total_heat_is_nonzero_and_balanced(air):
    hx = make_exchanger(air, counter_current=False)
    hx["hot_side"]["mass_flow"] = 2.0
    hx["cold_side"]["mass_flow"] = 2.0 / hx["mass_flow_ratio"]
    compute_component_energy_flows({"hx": hx})
    assert hx["heat_flow"] == pytest.approx(2.0 * hx["q_hot_side"])
    assert hx["heat_balance"] == pytest.approx(0.0, abs=1e-7)


@pytest.mark.parametrize("process,machine,p_in,t_in,p_out,temperature,work", [
    # Compressor reference includes the Figure 4 correction to printed Eq. 23.
    (compression_process, compressor_defaults, 1e5, 300, 2e5, 378.3336750060007, 78876.53049150406),
    (expansion_process, turbine_defaults, 2e5, 600, 1e5, 507.55686502716463, 96219.04336119222),
], ids=["centrifugal_compressor", "radial_turbine"])
def test_nondimensional_machine_outlet_and_design_regression(air, process, machine, p_in, t_in, p_out, temperature, work):
    inlet = air.get_state(props.PT_INPUTS, p_in, t_in)
    component = process(air, inlet.h, p_in, p_out, efficiency_type="non-dimensional",
                        mass_flow=2.0, data_in=machine())
    assert_finite_state(component["state_out"])
    assert component["state_out"].T == pytest.approx(temperature, rel=2e-6)
    assert component["specific_work"] == pytest.approx(work, rel=2e-6)
    assert component["state_out"].p == pytest.approx(p_out)
    output = component["data_out"]
    assert output["power"] == pytest.approx(2.0 * work, rel=2e-6)
    for field in ("diameter", "peripheral_speed", "angular_speed", "mechanical_stress"):
        assert np.isfinite(output[field]) and output[field] > 0
    assert 0 < output["isentropic_efficiency"] < 1
    matched_efficiency = output["polytropic_efficiency"] if process is compression_process else output["isentropic_efficiency"]
    assert output["correlation_efficiency"] == pytest.approx(matched_efficiency, abs=1e-10)


def test_shrouded_compressor_can_evaluate_an_outlet(air):
    data = compressor_defaults()
    data["is_shrouded"] = True
    inlet = air.get_state(props.PT_INPUTS, 1e5, 300.0)
    component = compression_process(air, inlet.h, inlet.p, 2e5, efficiency_type="non-dimensional",
                                    mass_flow=2.0, data_in=data)
    assert_finite_state(component["state_out"])


def test_friction_correlation_regression_and_roughness_trend():
    assert friction_factor(1e5, 1e-5) == pytest.approx(0.03219134303488542, rel=1e-10)
    assert friction_factor(1e6, 1e-5) < friction_factor(1e5, 1e-5)
    assert friction_factor(1e6, 1e-3) > friction_factor(1e6, 1e-5) > 0
