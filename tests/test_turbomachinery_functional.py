"""Numerical edge cases and JAX transformations of functional machine models."""

import jax
import jax.numpy as jnp
import jaxprop as props
import numpy as np
import pytest

from thermopt.components.turbomachinery_nondimensional import (
    compressor_defaults, compressor_performance, polytropic_head,
    solve_compressor, solve_turbine, turbine_performance,
)


@pytest.mark.parametrize("flow,coefficient", [(.02, .67), (.05, .67), (.085, .71), (.12, .75), (.14, .75)])
def test_shrouded_stress_interpolation_and_clipped_endpoints(flow, coefficient):
    """Appendix B's stress interpolation must work inside and outside its span."""
    fluid = props.Fluid("Air", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, 1e5, 300.)
    result = solve_compressor(fluid, inlet, 2e5, 2., dict(is_shrouded=True, flow_coefficient=flow))
    output = result["performance"]
    assert output["mechanical_stress"] == pytest.approx(coefficient * 7800 * output["peripheral_speed"]**2)
    assert 0 < output["polytropic_efficiency"] < 1
    assert abs(result["residual"]) < 1e-10


@pytest.mark.parametrize("offset", [0., 1e-10, -1e-10, .02, -.02])
def test_polytropic_head_is_finite_and_differentiable_at_equal_temperatures(offset):
    """The logarithmic-mean limit is T, with temperature derivative 1/2."""
    inlet = dict(temperature=300., enthalpy=4e5, entropy=2000.)

    def head(temperature):
        """Evaluate Eq. 9 with fixed h,s to isolate the temperature mean."""
        return polytropic_head(inlet, dict(temperature=temperature, enthalpy=4.1e5, entropy=2001.), xp=jnp)

    value, derivative = jax.jit(jax.value_and_grad(head))(300. + offset)
    change = offset / 300.
    expected_mean = 300. * (1 + change / 2 - change**2 / 12 + change**3 / 24)
    assert value == pytest.approx(10000. - expected_mean, abs=1e-9)
    assert derivative == pytest.approx(-.5 + change / 6 - change**2 / 8, abs=1e-8)


@pytest.mark.parametrize("compressor", [True, False], ids=["compressor", "turbine"])
def test_jit_grad_and_vmap_match_numpy_and_finite_differences(compressor):
    """Trace complete trial-state evaluation, including a perfect-gas h,s flash."""
    fluid = props.FluidPerfectGas("Air")
    p_in, temperature, p_out = (1e5, 300., 2e5) if compressor else (2e5, 600., 1e5)
    inlet = fluid.get_state(props.PT_INPUTS, p_in, temperature)
    isentropic = fluid.get_state(props.PSmass_INPUTS, p_out, inlet.s)
    work = isentropic.h - inlet.h
    enthalpy = inlet.h + work * (1 / .85 if compressor else .85)
    point = np.array([.075 if compressor else .55, float(enthalpy)])
    fields = ["correlation_efficiency", "polytropic_efficiency", "isentropic_efficiency",
              "diameter", "angular_speed", "mechanical_stress"]

    def evaluate(parameters, xp):
        """Evaluate geometry and efficiencies for a design coordinate and outlet h."""
        outlet = fluid.get_state(props.HmassP_INPUTS, parameters[1], p_out, errors="mask")
        if compressor:
            result = compressor_performance(inlet, outlet, isentropic, 2.,
                                            {"flow_coefficient": parameters[0]}, xp=xp)
        else:
            result = turbine_performance(fluid, inlet, outlet, isentropic, 2.,
                                         {"specific_speed": parameters[0]}, xp=xp, errors="mask")
        return xp.stack([result[field] for field in fields])

    def traced(parameters):
        """Select JAX arithmetic while keeping the model and fluid explicit."""
        return evaluate(parameters, jnp)

    np.testing.assert_allclose(jax.jit(traced)(point), evaluate(point, np), rtol=1e-11)
    jacobian = jax.jit(jax.jacfwd(traced))(point)
    finite_difference = []
    for column, step in enumerate([1e-6, .05]):
        delta = np.zeros(2)
        delta[column] = step
        finite_difference.append((evaluate(point + delta, np) - evaluate(point - delta, np)) / (2 * step))
    np.testing.assert_allclose(jacobian, np.array(finite_difference).T, rtol=2e-5, atol=1e-8)
    batch = np.stack([point, point + [.001, 50.], point - [.001, 50.]])
    np.testing.assert_allclose(jax.jit(jax.vmap(traced))(batch),
                               np.stack([evaluate(item, np) for item in batch]), rtol=1e-11)


@pytest.mark.parametrize("solve,p_in,temperature,p_out", [
    (solve_compressor, 1e5, 300., 2e5), (solve_turbine, 2e5, 600., 1e5)])
def test_host_solver_accepts_perfect_gas_fluid(solve, p_in, temperature, p_out):
    """A supplied perfect-gas model also satisfies the scalar solver contract."""
    fluid = props.FluidPerfectGas("Air")
    inlet = fluid.get_state(props.PT_INPUTS, p_in, temperature)
    result = solve(fluid, inlet, p_out, 2.)
    assert abs(result["residual"]) < 1e-10
    assert result["state_out"].p == pytest.approx(p_out)


@pytest.mark.parametrize("design,message", [
    ({"flow_coefficient": 0}, "flow_coefficient"), ({"roughness": -1}, "roughness"),
    ({"metal_density": 0}, "metal_density"), ({"backsweep": np.nan}, "backsweep"),
    ({"has_vaned": 1}, "has_vaned"), ({"unknown": 1}, "Unknown design"),
])
def test_bad_compressor_design_is_rejected(design, message):
    """Invalid inputs fail at the public boundary, before a numerical failure."""
    fluid = props.Fluid("Air", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, 1e5, 300.)
    with pytest.raises(ValueError, match=message):
        solve_compressor(fluid, inlet, 2e5, 2., design)


@pytest.mark.parametrize("mass,p_out,message", [
    (0, 2e5, "mass_flow"), (np.inf, 2e5, "mass_flow"),
    (None, 2e5, "mass_flow"), (2., np.nan, "p_out"),
    (2., 1e5, "greater"), (2., .9e5, "greater"),
])
def test_invalid_flow_and_compressor_pressure_direction(mass, p_out, message):
    """Degenerate and reversed compressor processes have no design solution."""
    fluid = props.Fluid("Air", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, 1e5, 300.)
    with pytest.raises(ValueError, match=message):
        solve_compressor(fluid, inlet, p_out, mass)


@pytest.mark.parametrize("design,p_out,message", [
    ({"specific_speed": 0}, 1e5, "specific_speed"),
    ({"diffuser_loss_fraction": 1.1}, 1e5, "diffuser_loss_fraction"),
    ({}, 3e5, "less"),
])
def test_invalid_turbine_design_and_pressure_direction(design, p_out, message):
    """Reject invalid turbine design domains and a pressure increase."""
    fluid = props.Fluid("Air", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, 2e5, 600.)
    with pytest.raises(ValueError, match=message):
        solve_turbine(fluid, inlet, p_out, 2., design)


def test_unbracketed_design_reports_context():
    """A deliberately impossible correction should identify the failing solve."""
    fluid = props.Fluid("Air", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, 1e5, 300.)
    with pytest.raises(ValueError, match="Could not solve compressor.*p_out=200000.*mass_flow=2"):
        solve_compressor(fluid, inlet, 2e5, 2., {"d_eta_other": -2.})


def test_default_dictionaries_do_not_share_state():
    """Independent designs must not inherit a previous caller's mutations."""
    original = compressor_defaults()
    changed = compressor_defaults()
    changed["roughness"] = 1.
    assert compressor_defaults() == original
