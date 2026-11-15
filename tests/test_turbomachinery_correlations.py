"""Independent paper checks and frozen radial-machine regression cases.

Paper: https://doi.org/10.1016/j.enconman.2025.120085
Table 4 / Eqs. 12-27 and Figure 6 / Eqs. 29-41 (Table 6 has a typo).
Compressor regression references include the approved Figure 4 correction
to the erroneous printed Equation 23. Equation checks do not constitute
validation against the authors' mean-line or experimental datasets.
"""

import json
from math import log, pi
from pathlib import Path

import jaxprop as props
import numpy as np
import pytest

from thermopt.components.turbomachinery_nondimensional import (
    compressor_performance, friction_factor, solve_compressor, solve_turbine,
)


@pytest.mark.parametrize("phi", [.02, .07, .14])
@pytest.mark.parametrize("backsweep", [0., 20., 40.])
@pytest.mark.parametrize("vaned", [False, True])
def test_compressor_equations_at_reference_size_speed_and_clearance(phi, backsweep, vaned):
    # Set the thermodynamic input states so D=D_ref, U=U_ref, Ma=0.8.
    # Both Reynolds and Mach corrections then vanish independently of phi.
    work_coefficient = .68 + (phi / (-pi * .37 / 4)) ** 3 + pi * .002 / (4 * phi)
    head_coefficient = .59 + 4 * .7 * phi / pi - (4 / pi) ** 2 * 7.5 * phi**2 - pi * .00025 / (4 * phi)
    efficiency = head_coefficient / work_coefficient
    if not vaned:
        efficiency -= .017 / (.04 + 20 * phi / pi + efficiency**3)
    efficiency -= .0006 * (40 - backsweep)
    actual_work_coefficient = work_coefficient + (.004 + .2 * phi**2) * (40 - backsweep)
    work = actual_work_coefficient * 277**2
    design = dict(flow_coefficient=phi, backsweep=backsweep, has_vaned=vaned)
    rho = 1.2
    mass_flow = rho * phi * 277 * .4**2
    inlet = dict(enthalpy=4e5, entropy=2000., temperature=300., density=rho,
                         viscosity=rho * 15.7e-6, speed_of_sound=277 / .8)
    log_mean_temperature = 100 / log(400 / 300)
    outlet = dict(enthalpy=4e5 + work, entropy=2000 + .2 * work / log_mean_temperature, temperature=400.)
    outlet_isentropic = dict(enthalpy=4e5 + .75 * work)
    output = compressor_performance(inlet, outlet, outlet_isentropic, mass_flow, design)
    residual = output["correlation_efficiency"] - output["polytropic_efficiency"]
    assert output["correlation_efficiency"] == pytest.approx(efficiency, abs=1e-12)
    assert residual == pytest.approx(efficiency - .8, abs=1e-12)
    assert output["diameter"] == pytest.approx(.4)
    assert output["peripheral_speed"] == pytest.approx(277)
    assert output["mechanical_stress"] == pytest.approx(.41 * 7800 * 277**2)


@pytest.mark.parametrize("phi,velocity_ratio", [(.01, .40922), (.07, .55178), (.15, .6545)])
@pytest.mark.parametrize("diameter,tip_speed", [(.04, 200.), (.4, 500.)])
def test_reynolds_correction_uses_figure4_velocity_for_actual_and_reference_states(
    phi, velocity_ratio, diameter, tip_speed,
):
    """Check the corrected velocity's coupling to both chord Reynolds numbers.

    The velocity ratios are independent Figure 4 values. Construct a known
    size/speed with zero Mach and clearance penalties, so an incorrect actual
    or reference Reynolds velocity cannot be hidden by other corrections.
    """
    work_coefficient = .68 + (phi / (-pi * .37 / 4))**3 + pi * .002 / (4 * phi)
    head_coefficient = .59 + 4 * .7 * phi / pi - (4 / pi)**2 * 7.5 * phi**2 - pi * .00025 / (4 * phi)
    work = work_coefficient * tip_speed**2
    rho, viscosity = 1.2, 15.7e-6
    inlet = dict(enthalpy=4e5, entropy=2000., temperature=300., density=rho, viscosity=rho * viscosity, speed_of_sound=tip_speed / .8)
    outlet = dict(enthalpy=4e5 + work, entropy=2000. + .2 * work / (100 / log(400 / 300)), temperature=400.)
    isentropic = dict(enthalpy=4e5 + .75 * work)
    mass_flow = rho * phi * tip_speed * diameter**2
    output = compressor_performance(inlet, outlet, isentropic, mass_flow, {"flow_coefficient": phi})
    diameter_over_chord = 1.1 + 22 * phi
    chord = diameter / diameter_over_chord
    reference_chord = .4 / diameter_over_chord
    actual_reynolds = velocity_ratio * tip_speed * chord / viscosity
    reference_reynolds = velocity_ratio * 277 * reference_chord / 15.7e-6
    actual_friction = friction_factor(actual_reynolds, 3.6e-6 / chord)
    reference_friction = friction_factor(reference_reynolds, 3.6e-6 / reference_chord)
    expected_correction = -(.05 + .002 / (phi + .0025)) * (actual_friction - reference_friction) / reference_friction
    assert output["correlation_efficiency"] == pytest.approx(
        head_coefficient / work_coefficient + expected_correction, abs=1e-12,
    )


@pytest.mark.parametrize("speed", [.35, .45, .60, .85])
def test_turbine_efficiency_agrees_with_published_figure_6(speed):
    fluid = props.Fluid("Air", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, 2e5, 600.)
    result = solve_turbine(fluid, inlet, 1e5, 2., {"specific_speed": speed})
    # Figure 6 and the author code take precedence over Table 6's 0.30 typo.
    base = .93 - .30 * (speed - .6)**2 - .50 * (speed - .6)**3
    if speed < .45:
        base -= 4.2 * (speed - .45)**2
    rotor_efficiency = base - .03
    kinetic_fraction = max(.02 + .08 * speed**2, .2 * speed**4)
    overall = rotor_efficiency / (1 + .35 * kinetic_fraction)
    output = result["performance"]
    reynolds = 4 / (inlet.mu * output["diameter"])
    expected = 1 - (1 - overall) * (.3 + .7 * min(1, reynolds / 352000.)**(-.2))
    assert output["correlation_efficiency"] == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("case", json.loads(
    (Path(__file__).parent / "data/turbomachinery_real_fluid.json").read_text()))
def test_enthalpy_solver_matches_high_pressure_air_and_co2_regression(case):
    """Check four frozen designs; compressor references include the Figure 4 fix."""
    fluid = props.Fluid(case["fluid"], backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, case["p_in"], case["T_in"])
    solve = solve_compressor if case["kind"] == "compressor" else solve_turbine
    result = solve(fluid, inlet, case["p_out"], case["mass_flow"])
    aliases = dict(T="T", p="p", hmass="h", smass="s", rhomass="rho",
                   viscosity="mu", speed_sound="a")
    for key, expected in case["outlet"].items():
        assert result["state_out"][aliases[key]] == pytest.approx(expected, rel=1e-8, abs=1e-7)
    for key, expected in case["performance"].items():
        assert result["performance"][key] == pytest.approx(expected, rel=1e-8, abs=1e-7)
    assert abs(result["residual"]) < 1e-10


@pytest.fixture(scope="module")
def demo_rows():
    from development.turbomachinery import demo_compressor, demo_turbine
    return {"compressor": demo_compressor.run_sweep(points=5),
            "turbine": demo_turbine.run_sweep(points=5)}


@pytest.mark.parametrize("kind", ["compressor", "turbine"])
def test_development_sweeps_match_recorded_outputs(kind, demo_rows):
    """Protect recorded sweeps after the deliberate compressor model correction."""
    reference = json.loads((Path(__file__).parent / "data/turbomachinery_trends.json").read_text())
    rows = demo_rows[kind]
    expected = reference[kind]
    assert len(rows) == len(expected)
    for actual, saved in zip(rows, expected):
        assert actual.keys() == saved.keys()
        for key, value in saved.items():
            if isinstance(value, str):
                assert actual[key] == value
            else:
                assert np.isfinite(actual[key])
                assert actual[key] == pytest.approx(value, rel=2e-7, abs=1e-7), key
        assert abs(actual["residual"]) < 1e-10
        for key in ["polytropic_efficiency", "isentropic_efficiency"]:
            assert 0 < actual[key] < 1
        for key in ["diameter", "power", "mechanical_stress", "angular_speed"]:
            assert actual[key] > 0


def curve(rows, family, label, output):
    """Extract one ordered output curve from the flat demo rows."""
    return np.array([row[output] for row in rows if row["family"] == family and row["case"] == label])


def test_compressor_trends_have_expected_design_and_scale_effects(demo_rows):
    rows = demo_rows["compressor"]
    output = "polytropic_efficiency"
    baseline = curve(rows, "design choices", "vaned, backsweep 40 deg", output)
    for label in ["vaneless, backsweep 40 deg", "vaned, backsweep 0 deg", "vaned, clearance 4%"]:
        assert np.all(curve(rows, "design choices", label, output) < baseline)
    small, medium, large = [curve(rows, "mass flow", f"m={mass:g} kg/s", output) for mass in [.1, 2., 20.]]
    assert np.all(small < medium)
    assert np.all(medium < large)
    peaks = [max(curve(rows, "pressure ratio", f"PR={pr}", output)) for pr in [2, 4, 6]]
    assert peaks[0] > peaks[1] > peaks[2]


def test_turbine_losses_and_scale_have_expected_effects(demo_rows):
    rows = demo_rows["turbine"]
    output = "isentropic_efficiency"
    baseline = curve(rows, "design choices", "scalloped, diffuser loss 0.35", output)
    assert np.all(curve(rows, "design choices", "unscalloped, diffuser loss 0.35", output) > baseline)
    for label in ["scalloped, diffuser loss 0.50", "scalloped, clearance 4%"]:
        assert np.all(curve(rows, "design choices", label, output) < baseline)
    assert np.all(curve(rows, "mass flow", "m=2 kg/s", output) >= curve(rows, "mass flow", "m=0.1 kg/s", output))


def test_turbine_scale_invariance_above_reynolds_threshold():
    fluid = props.Fluid("Air", backend="HEOS")
    inlet = fluid.get_state(props.PT_INPUTS, 2e5, 600.)
    small, large = [solve_turbine(fluid, inlet, 1e5, mass)["performance"] for mass in [2., 20.]]
    assert large["isentropic_efficiency"] == pytest.approx(small["isentropic_efficiency"], abs=1e-10)
    assert large["diameter"] / small["diameter"] == pytest.approx(np.sqrt(10))
    assert large["angular_speed"] / small["angular_speed"] == pytest.approx(1 / np.sqrt(10))
    assert large["power"] / small["power"] == pytest.approx(10)
