"""Polytropic path boundaries fail clearly and preserve complete reporting."""

import jaxprop as jp
import numpy as np
import pytest

from thermopt.components import basic_components as basic


@pytest.fixture(scope="module", params=["coolprop", "perfect_gas", "bicubic"])
def model(request, test_workspace):
    settings = dict(name="Air", method=request.param, identifier="boundary_air")
    if request.param == "bicubic":
        settings["table_options"] = dict(
            h_min=300e3,
            h_max=850e3,
            p_min=40e3,
            p_max=350e3,
            N_h=16,
            N_p=16,
            N_p_sat=4,
            coarse_step=1,
            table_dir=str(test_workspace / "boundary_tables"),
        )
    return jp.Fluid(**settings)


@pytest.mark.parametrize(
    "process", [basic.compression_process, basic.expansion_process]
)
@pytest.mark.parametrize("steps", [0, 1, -1, 2.0, np.nan, True])
def test_invalid_polytropic_output_count_fails_before_flashing(
    process, steps, monkeypatch
):
    model = jp.Fluid("Air", method="perfect_gas")

    def unexpected(*args, **kwargs):
        raise AssertionError("Invalid path count must fail at the component boundary")

    monkeypatch.setattr(type(model), "get_state", unexpected)
    with pytest.raises(ValueError, match="num_steps.*integer.*2"):
        process(
            model,
            450e3,
            1e5,
            2e5,
            efficiency=0.85,
            efficiency_type="polytropic",
            num_steps=steps,
        )


@pytest.mark.parametrize(
    "process", [basic.compression_process, basic.expansion_process]
)
@pytest.mark.parametrize("steps", [np.int64(2), 9])
def test_zero_pressure_span_returns_constant_full_state_path(
    model, process, steps, monkeypatch
):
    inlet = model.get_state(jp.PT_INPUTS, 1e5, 400.0)

    def unexpected(*args, **kwargs):
        raise AssertionError("Zero pressure span must not invoke integration")

    monkeypatch.setattr(basic, "solve_ivp", unexpected)
    result = process(
        model,
        float(inlet.enthalpy),
        1e5,
        1e5,
        efficiency=0.85,
        efficiency_type="polytropic",
        num_steps=steps,
    )
    path = result["states"]
    assert path.shape == (steps,)
    assert (
        path.identifier == model.identifier and path.reference_id == model.reference_id
    )
    np.testing.assert_array_equal(path.status, 0)
    np.testing.assert_allclose(path.enthalpy, inlet.enthalpy, atol=1e-5)
    np.testing.assert_allclose(path.pressure, 1e5)
    np.testing.assert_allclose(path.temperature, 400.0, atol=1e-6)
    assert result["specific_work"] == pytest.approx(0.0, abs=1e-8)
    assert result["isentropic_work"] == pytest.approx(0.0, abs=1e-5)
    for end in ("state_in", "state_out"):
        assert result[end].shape == ()
        assert result[end].enthalpy == pytest.approx(inlet.enthalpy, abs=1e-5)
    if model.provider == "coolprop":
        assert "superheating" in path.available_properties


@pytest.mark.parametrize(
    "process", [basic.compression_process, basic.expansion_process]
)
def test_unreturned_isentropic_reference_omits_optional_diagnostics(
    process, monkeypatch
):
    model = jp.Fluid("CO2")
    p_in, p_out, temp = (
        (1e6, 2e6, 300.0) if process is basic.compression_process else (2e6, 1e6, 600.0)
    )
    inlet = model.get_state(jp.PT_INPUTS, p_in, temp)
    calls = []
    original = type(model).get_state

    def observed(self, pair, first, second, **options):
        calls.append((pair, options))
        return original(self, pair, first, second, **options)

    monkeypatch.setattr(type(model), "get_state", observed)
    result = process(model, float(inlet.enthalpy), p_in, p_out, efficiency=0.85)
    assert [options for pair, options in calls if pair == jp.PSmass_INPUTS] == [{}]
    assert np.isfinite(result["state_in"].superheating)
    assert np.isfinite(result["state_out"].superheating)


@pytest.mark.parametrize(
    "process", [basic.compression_process, basic.expansion_process]
)
@pytest.mark.parametrize("value", [True, "450000", 450e3 + 1j])
def test_host_boundary_preserves_native_invalid_input_rejection(process, value):
    model = jp.Fluid("Air", method="perfect_gas")
    with pytest.raises(jp.InvalidStateError):
        process(model, value, 1e5, 2e5, efficiency=0.85, efficiency_type="polytropic")


@pytest.mark.parametrize(
    "process", [basic.compression_process, basic.expansion_process]
)
def test_prepared_state_shapes_do_not_recompile_for_scipy_scalar_types(process, caplog):
    """Guard compilation reuse without asserting shared-machine wall times."""
    import jax

    compression = process is basic.compression_process
    model = jp.Fluid("Air", method="perfect_gas", identifier="polytropic_scalar_cache")
    p_in, temperature, p_out = (1e5, 300.0, 2e5) if compression else (2e5, 600.0, 1e5)
    inlet = model.get_state(jp.PT_INPUTS, p_in, temperature)
    h = float(inlet.enthalpy)
    scalar = model.get_state(jp.HmassP_INPUTS, np.asarray(h), np.asarray(p_in))
    isentropic = model.get_state(jp.PSmass_INPUTS, np.asarray(p_out), scalar.entropy)
    path = model.get_state(jp.HmassP_INPUTS, np.full(9, h), np.linspace(p_in, p_out, 9))
    jax.block_until_ready((scalar, isentropic, path))
    caplog.clear()
    with jax.log_compiles(True):
        result = process(
            model,
            h,
            p_in,
            p_out,
            efficiency=0.85,
            efficiency_type="polytropic",
            num_steps=9,
        )
        jax.block_until_ready(result)
    compilations = [
        record.getMessage()
        for record in caplog.records
        if "Compiling get_state with" in record.getMessage()
    ]
    assert (
        not compilations
    ), "SciPy scalar representation caused a redundant state compilation"
    np.testing.assert_array_equal(result["states"].status, 0)
    assert result["state_out"].pressure == pytest.approx(p_out)
