"""Independent eager component physics through the native JaxProp models.

Perfect-gas references follow the constant-cp ideal-gas equations. Bicubic
round trips are checked inside its declared Air domain; they do not assert
universal phase/saturation interchangeability. Existing source-fluid numerical
fixtures retain their original tolerances in test_components.py.
"""

from copy import deepcopy
from types import SimpleNamespace

import jaxprop as jp
import numpy as np
import pytest

import thermopt as th
from thermopt.components import basic_components as basic
from thermopt.components.turbomachinery_nondimensional import (
    compressor_defaults,
    solve_compressor,
    solve_turbine,
    turbine_defaults,
)
from thermopt.utilities import process_state_options
from regression_support import load_case


@pytest.fixture(scope="module")
def models(test_workspace):
    settings = dict(name="Air", identifier="component_air")
    tables = dict(
        h_min=200e3,
        h_max=850e3,
        p_min=40e3,
        p_max=350e3,
        N_h=32,
        N_p=32,
        N_p_sat=4,
        coarse_step=1,
        table_dir=str(test_workspace / "component_tables"),
    )
    return {
        method: jp.Fluid(
            **settings,
            method=method,
            **(dict(table_options=tables) if method == "bicubic" else {})
        )
        for method in ("coolprop", "perfect_gas", "bicubic")
    }


@pytest.fixture(scope="module", params=["coolprop", "perfect_gas", "bicubic"])
def model(request, models):
    return models[request.param]


MACHINES = [
    (basic.compression_process, 1e5, 300.0, 2e5),
    (basic.expansion_process, 2e5, 600.0, 1e5),
]


def assert_physical_path(path, model):
    assert path.identifier == model.identifier
    assert path.reference_id == model.reference_id
    assert path.provider == model.provider
    np.testing.assert_array_equal(path.status, 0)
    for field in ("temperature", "pressure", "density", "cp", "cv", "speed_sound"):
        assert np.all(np.isfinite(path[field])), field
        assert np.all(np.asarray(path[field]) > 0), field
    for field in ("enthalpy", "entropy"):
        assert np.all(np.isfinite(path[field])), field


@pytest.mark.parametrize("process,p_in,t_in,p_out", MACHINES)
@pytest.mark.parametrize("kind", ["isentropic", "polytropic", "non-dimensional"])
def test_machinery_energy_efficiency_and_inputs(
    model, process, p_in, t_in, p_out, kind
):
    inlet = model.get_state(jp.PT_INPUTS, p_in, t_in)
    compressor = process is basic.compression_process
    design = compressor_defaults() if compressor else turbine_defaults()
    original = deepcopy(design)
    result = process(
        model,
        float(inlet.enthalpy),
        p_in,
        p_out,
        efficiency=0.85,
        efficiency_type=kind,
        mass_flow=2.0,
        data_in=design,
        num_steps=9,
    )
    assert design == original
    assert_physical_path(result["states"], model)
    outlet = result["state_out"]
    assert outlet.pressure == pytest.approx(p_out, abs=1e-5)
    assert result["state_in"].temperature == pytest.approx(t_in, abs=0.02)
    sign = 1 if compressor else -1
    assert result["specific_work"] == pytest.approx(
        sign * (outlet.enthalpy - inlet.enthalpy), abs=1e-5
    )
    assert result["specific_work"] > 0
    assert outlet.entropy > inlet.entropy
    isentropic = model.get_state(jp.PSmass_INPUTS, p_out, inlet.entropy)
    eta_is = (isentropic.enthalpy - inlet.enthalpy) / (outlet.enthalpy - inlet.enthalpy)
    if not compressor:
        eta_is = 1 / eta_is
    assert 0 < eta_is < 1
    if kind == "isentropic":
        assert eta_is == pytest.approx(0.85, abs=1e-8)
    elif kind == "non-dimensional":
        solve = solve_compressor if compressor else solve_turbine
        solved = solve(model, inlet, p_out, 2.0, design)
        assert abs(solved["residual"]) < 1e-10
        assert result["data_out"]["power"] == pytest.approx(2 * result["specific_work"])
        assert result["data_out"]["isentropic_efficiency"] == pytest.approx(eta_is)
        assert outlet.enthalpy == pytest.approx(solved["state_out"].enthalpy)


@pytest.mark.parametrize("process,p_in,t_in,p_out", MACHINES)
def test_reversible_isentropic_limit_for_native_models(
    model, process, p_in, t_in, p_out
):
    inlet = model.get_state(jp.PT_INPUTS, p_in, t_in)
    result = process(model, float(inlet.enthalpy), p_in, p_out, efficiency=1.0)
    assert result["state_out"].entropy == pytest.approx(inlet.entropy, abs=0.002)
    assert result["specific_work"] == pytest.approx(result["isentropic_work"])


@pytest.mark.parametrize("counter_current", [False, True])
@pytest.mark.parametrize("zero_duty", [False, True])
def test_exchanger_heat_balance_and_valve_throttling(model, counter_current, zero_duty):
    hot_in, hot_out = [model.get_state(jp.PT_INPUTS, 2e5, t) for t in (650.0, 400.0)]
    cold_in, cold_out = [model.get_state(jp.PT_INPUTS, 1e5, t) for t in (300.0, 350.0)]
    hx = basic.heat_exchanger(
        model,
        hot_in.enthalpy,
        hot_in.enthalpy if zero_duty else hot_out.enthalpy,
        2e5,
        2e5,
        model,
        cold_in.enthalpy,
        cold_in.enthalpy if zero_duty else cold_out.enthalpy,
        1e5,
        1e5,
        num_steps=9,
        counter_current=counter_current,
    )
    for side in ("hot_side", "cold_side"):
        assert_physical_path(hx[side]["states"], model)
    # m_hot/m_cold = q_cold/q_hot, independently of plotting direction.
    if zero_duty:
        assert hx["q_hot_side"] == hx["q_cold_side"] == 0
        assert hx["mass_flow_ratio"] == 1
    else:
        assert hx["mass_flow_ratio"] * hx["q_hot_side"] == pytest.approx(
            hx["q_cold_side"]
        )
        hot_temperatures = [400.0, 650.0] if counter_current else [650.0, 400.0]
        np.testing.assert_allclose(
            np.asarray(hx["hot_side"]["states"].temperature)[[0, -1]],
            hot_temperatures,
            atol=0.02,
        )
    assert np.all(hx["temperature_difference"] > 0)
    valve = basic.isenthalpic_valve(hot_in, 1e5, model, N=7)
    assert_physical_path(valve, model)
    np.testing.assert_allclose(valve.enthalpy, hot_in.enthalpy, atol=1e-5)
    np.testing.assert_allclose(valve.pressure, np.linspace(2e5, 1e5, 7))
    assert np.all(np.diff(valve.entropy) > 0)


@pytest.mark.parametrize(
    "pair", [jp.HmassP_INPUTS, jp.PT_INPUTS, jp.PSmass_INPUTS, jp.HmassSmass_INPUTS]
)
def test_native_pair_round_trips(model, pair):
    target = model.get_state(jp.HmassP_INPUTS, 600e3, 150e3)
    values = {
        jp.HmassP_INPUTS: (target.enthalpy, target.pressure),
        jp.PT_INPUTS: (target.pressure, target.temperature),
        jp.PSmass_INPUTS: (target.pressure, target.entropy),
        jp.HmassSmass_INPUTS: (target.enthalpy, target.entropy),
    }[pair]
    recovered = model.get_state(pair, *values)
    assert_physical_path(recovered, model)
    np.testing.assert_allclose(recovered.temperature, target.temperature, atol=0.002)
    np.testing.assert_allclose(recovered.pressure, target.pressure, rtol=2e-5)


@pytest.mark.parametrize("h,p_in,p_out", [(300e3, 350e3, 40e3), (850e3, 350e3, 40e3)])
def test_bicubic_valve_includes_domain_edges_and_rejects_outside(
    models, h, p_in, p_out
):
    model = models["bicubic"]
    inlet = model.get_state(jp.HmassP_INPUTS, h, p_in)
    path = basic.isenthalpic_valve(inlet, p_out, model, N=9)
    assert_physical_path(path, model)
    np.testing.assert_allclose(path.enthalpy, h)
    np.testing.assert_allclose(np.asarray(path.pressure)[[0, -1]], [p_in, p_out])
    with pytest.raises(jp.OutOfDomainError):
        basic.isenthalpic_valve(inlet, p_out - 100.0, model, N=9)


def test_bicubic_low_enthalpy_corner_distinguishes_valid_gas_from_wet_cell(models):
    table, source = models["bicubic"], models["coolprop"]
    # The same declared h_min is a supported gas state at p_min and a wet
    # Air state at p_max. A rectangular bounding box is not a trust guarantee.
    assert_physical_path(table.get_state(jp.HmassP_INPUTS, 200e3, 40e3), table)
    assert source.get_state(jp.HmassP_INPUTS, 200e3, 350e3).is_two_phase == 1
    with pytest.raises(jp.OutOfDomainError, match="trusted interpolation cell"):
        table.get_state(jp.HmassP_INPUTS, 200e3, 350e3)


@pytest.mark.parametrize("method", ["perfect_gas", "bicubic"])
def test_saturation_is_explicitly_model_specific(models, method):
    model = models[method]
    with pytest.raises(jp.UnsupportedInputPairError):
        model.get_state(jp.PQ_INPUTS, 1e5, 0.5)
    with pytest.raises(jp.UnsupportedOperationError, match="saturation"):
        basic.heat_exchanger(
            model,
            700e3,
            500e3,
            2e5,
            2e5,
            model,
            300e3,
            400e3,
            1e5,
            1e5,
            include_saturation_nodes=True,
        )


@pytest.mark.parametrize("process,p_in,t_in,p_out", MACHINES)
@pytest.mark.parametrize("kind", ["isentropic", "polytropic"])
def test_perfect_gas_machinery_matches_analytic_temperature(
    process, p_in, t_in, p_out, kind
):
    model = jp.Fluid("Air", method="perfect_gas")
    inlet = model.get_state(jp.PT_INPUTS, p_in, t_in)
    gamma, cp = float(inlet.cp / inlet.cv), float(inlet.cp)
    exponent = (gamma - 1) / gamma
    compression = process is basic.compression_process
    eta = 0.85
    if kind == "polytropic":
        expected = t_in * (p_out / p_in) ** (
            exponent / eta if compression else exponent * eta
        )
    else:
        reversible = t_in * (p_out / p_in) ** exponent
        expected = t_in + (reversible - t_in) * (1 / eta if compression else eta)
    result = process(
        model,
        float(inlet.enthalpy),
        p_in,
        p_out,
        efficiency=eta,
        efficiency_type=kind,
        num_steps=9,
    )
    # RK45 defaults have rtol=1e-3. Preserve those solver settings; the analytic
    # check is separate from the much tighter frozen source regression checks.
    tolerance = 1e-4 if kind == "polytropic" else 1e-11
    assert result["state_out"].temperature == pytest.approx(expected, rel=tolerance)
    assert result["specific_work"] == pytest.approx(
        cp * abs(expected - t_in), rel=5 * tolerance
    )


@pytest.mark.parametrize("kind", ["isentropic", "polytropic"])
def test_complete_eager_perfect_gas_cycle_obeys_ideal_gas_first_law(kind, tmp_path):
    config = deepcopy(load_case("sco2_simple_polytropic")["problem_formulation"])
    params = config["fixed_parameters"]
    params["working_fluid"] = dict(
        name="Air", method="perfect_gas", model_options=dict(T_ref=600.0, p_ref=1e5)
    )
    for machine in ("compressor", "expander"):
        params[machine]["efficiency_type"] = kind
    model = jp.Fluid(**params["working_fluid"])
    inlet_c = model.get_state(jp.PT_INPUTS, 1e5, 300.0)
    inlet_t = model.get_state(jp.PT_INPUTS, 2e5, 750.0)
    replacements = dict(
        compressor_inlet_pressure=1e5,
        expander_inlet_pressure=2e5,
        compressor_inlet_enthalpy=float(inlet_c.enthalpy),
        expander_inlet_enthalpy=float(inlet_t.enthalpy),
    )
    for name, value in replacements.items():
        config["design_variables"][name] = dict(
            min=0.8 * value, value=value, max=1.2 * value
        )
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    identities = {role: id(fluid) for role, fluid in problem.fluids.items()}
    first = np.asarray(problem.fitness(problem.x0))
    np.testing.assert_array_equal(first, problem.fitness(problem.x0))
    assert np.all(np.isfinite(first))
    assert identities == {role: id(fluid) for role, fluid in problem.fluids.items()}
    components = problem.cycle_data["components"]
    c, t = components["compressor"], components["expander"]
    cp = float(inlet_c.cp)
    expected_c = cp * (c["state_out"].temperature - c["state_in"].temperature)
    expected_t = cp * (t["state_in"].temperature - t["state_out"].temperature)
    assert c["specific_work"] == pytest.approx(expected_c)
    assert t["specific_work"] == pytest.approx(expected_t)
    heater = components["heater"]["cold_side"]
    q_in = cp * (heater["state_out"].temperature - heater["state_in"].temperature)
    cooler = components["cooler"]["hot_side"]
    q_out = cp * (cooler["state_in"].temperature - cooler["state_out"].temperature)
    assert q_in - q_out == pytest.approx(expected_t - expected_c)
    assert problem.cycle_data["energy_analysis"]["energy_balance"] == pytest.approx(
        0.0, abs=1e-7
    )
    assert problem.cycle_data["energy_analysis"]["net_cycle_power"] == pytest.approx(
        1e6
    )


@pytest.mark.parametrize(
    "process,p_in,t_in,p_out",
    [
        (basic.compression_process, 1e6, 300.0, 2e6),
        (basic.expansion_process, 2e6, 600.0, 1e6),
    ],
)
def test_polytropic_reporting_matches_full_scalar_states(process, p_in, t_in, p_out):
    model = jp.Fluid("CO2", identifier="reported_source")
    inlet = model.get_state(jp.PT_INPUTS, p_in, t_in)
    result = process(
        model,
        float(inlet.enthalpy),
        p_in,
        p_out,
        efficiency=0.85,
        efficiency_type="polytropic",
        num_steps=9,
    )
    path = result["states"]
    # Pressures lie between CO2's triple and critical pressures, so saturation
    # reporting includes finite values as well as legitimately unavailable fields.
    assert np.all(np.isfinite(path.superheating))
    options = process_state_options(
        model, generalize_quality=process is basic.expansion_process
    )
    scalars = jp.FluidState.stack(
        [
            model.get_state(jp.HmassP_INPUTS, float(h), float(p), **options)
            for h, p in zip(path.enthalpy, path.pressure)
        ]
    )
    assert path.available_properties == scalars.available_properties
    actual, expected = path.to_data(), scalars.to_data()
    assert actual["metadata"] == expected["metadata"]
    for field in path.available_properties:
        np.testing.assert_allclose(
            actual["properties"][field],
            expected["properties"][field],
            rtol=1e-12,
            atol=1e-10,
            equal_nan=True,
            err_msg=field,
        )
        np.testing.assert_array_equal(
            actual["property_validity"][field], expected["property_validity"][field]
        )
    np.testing.assert_array_equal(path.status, scalars.status)
    assert path.reference_id == scalars.reference_id


@pytest.mark.parametrize("process,p_in,t_in,p_out", MACHINES)
def test_failed_polytropic_integration_is_not_returned_as_a_complete_component(
    monkeypatch, process, p_in, t_in, p_out
):
    model = jp.Fluid("Air", method="perfect_gas")
    inlet = model.get_state(jp.PT_INPUTS, p_in, t_in)
    monkeypatch.setattr(
        basic,
        "solve_ivp",
        lambda *args, **kwargs: SimpleNamespace(
            success=False,
            message="Required step size is less than spacing between numbers.",
            t=np.array([p_in, 0.5 * (p_in + p_out)]),
            y=np.array([[float(inlet.enthalpy), float(inlet.enthalpy)]]),
        ),
    )
    with pytest.raises(RuntimeError, match="[Pp]olytropic.*Required step size"):
        process(
            model,
            float(inlet.enthalpy),
            p_in,
            p_out,
            efficiency=0.85,
            efficiency_type="polytropic",
        )


@pytest.mark.parametrize("process,p_in,t_in,p_out", MACHINES)
def test_polytropic_trials_defer_reporting_to_complete_batched_path(
    monkeypatch, process, p_in, t_in, p_out
):
    model = jp.Fluid("CO2", identifier="reported_source")
    inlet = model.get_state(jp.PT_INPUTS, p_in, t_in)
    calls = []
    original = type(model).get_state

    def recorded(self, pair, first, second, **options):
        if self is model:
            calls.append((pair, np.shape(first), dict(options)))
        return original(self, pair, first, second, **options)

    monkeypatch.setattr(type(model), "get_state", recorded)
    result = process(
        model,
        float(inlet.enthalpy),
        p_in,
        p_out,
        efficiency=0.85,
        efficiency_type="polytropic",
        num_steps=9,
    )
    hp_calls = [
        (shape, options) for pair, shape, options in calls if pair == jp.HmassP_INPUTS
    ]
    reporting = process_state_options(
        model, generalize_quality=process is basic.expansion_process
    )
    assert hp_calls[0] == ((), reporting)  # Public inlet reporting is unchanged.
    assert hp_calls[-1] == ((9,), reporting)
    assert len(hp_calls[1:-1]) > 0
    assert all(shape == () and not options for shape, options in hp_calls[1:-1])
    assert "superheating" in result["states"].available_properties
