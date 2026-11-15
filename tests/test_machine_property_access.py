"""Verify supplied-fluid access, enthalpy-only iteration and rotor-static flashes."""

import jaxprop as props
import numpy as np
import pytest

from thermopt.components.turbomachinery_nondimensional import solve_compressor, solve_turbine


class RecordingFluid:
    """Record public property requests while delegating to a jaxprop model."""

    def __init__(self):
        """Create an independent Air model and an empty request log."""
        self.fluid = props.Fluid("Air", backend="HEOS")
        self.calls = []

    def get_state(self, pair, first, second):
        """Record each pair and its returned state without exposing the backend."""
        state = self.fluid.get_state(pair, first, second)
        self.calls.append((pair, first, second, state))
        return state


@pytest.mark.parametrize("solve,p_in,t_in,p_out", [
    (solve_compressor, 1e5, 300, 2e5), (solve_turbine, 2e5, 600, 1e5)])
def test_outlet_solver_reuses_inlet_and_only_flashes_enthalpy(solve, p_in, t_in, p_out):
    """The inlet is explicit; root evaluations use the same supplied fluid."""
    fluid = RecordingFluid()
    inlet = fluid.fluid.get_state(props.PT_INPUTS, p_in, t_in)
    before = {key: float(inlet[key]) for key in ["h", "s", "p", "T"]}
    design = {"metal_density": 8000.}
    result = solve(fluid, inlet, p_out, 2., design)
    pairs = [call[0] for call in fluid.calls]
    assert result["state_in"] is inlet
    assert pairs.count(props.PSmass_INPUTS) == 1
    assert props.PT_INPUTS not in pairs
    assert props.HmassP_INPUTS in pairs
    assert set(pairs) <= {props.PSmass_INPUTS, props.HmassP_INPUTS, props.HmassSmass_INPUTS}
    assert result["solver"]["converged"]
    assert result["solver"]["function_calls"] >= 3
    assert abs(result["residual"]) < 1e-10
    assert before == {key: float(inlet[key]) for key in before}
    assert design == {"metal_density": 8000.}


def test_rotor_static_outlet_uses_supplied_model_and_paper_entropy_allocation():
    """Check Eqs. 46–48's h,s flash independently at the converged outlet."""
    fluid = RecordingFluid()
    inlet = fluid.fluid.get_state(props.PT_INPUTS, 2e5, 600.)
    result = solve_turbine(fluid, inlet, 1e5, 2.)
    rotor_calls = [call for call in fluid.calls if call[0] == props.HmassSmass_INPUTS]
    _, enthalpy, entropy, rotor = rotor_calls[-1]
    outlet, isentropic = result["state_out"], result["state_out_isentropic"]
    kinetic_fraction = .02 + .08 * .6**2
    rotor_efficiency = .93 - .03
    overall = rotor_efficiency / (1 + .35 * kinetic_fraction)
    rotor_work = (isentropic.h - inlet.h) / (1 + .35 * kinetic_fraction)
    expected_entropy = inlet.s + (outlet.s - inlet.s) * (1 / rotor_efficiency - 1) / (1 / overall - 1)
    assert enthalpy == pytest.approx(outlet.h - kinetic_fraction * abs(rotor_work))
    assert entropy == pytest.approx(expected_entropy)
    assert np.isfinite(rotor.rho) and rotor.rho > 0
    expected_speed = .6 * abs(rotor_work)**.75 / np.sqrt(2 / rotor.rho)
    assert result["performance"]["angular_speed"] == pytest.approx(expected_speed)
