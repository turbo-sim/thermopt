"""Off-design tests. The end-to-end example runs only with THERMOPT_RUN_SLOW=1."""

import os
import runpy
import shutil
from pathlib import Path

import jaxprop as jxp
import numpy as np
import pytest

from thermopt.components import components_offdesign as cod

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "examples" / "rankine_partial_evaporation_offdesign"

# Turbine-inlet stagnation states for the nozzle tests: (fluid, pressure, quality)
NOZZLE_CASES = [("butane", 15e5, 0.05), ("butane", 15e5, 0.5), ("Novec649", 8e5, 0.05), ("Novec649", 8e5, 0.9)]


@pytest.fixture(scope="module")
def butane():
    return jxp.Fluid(name="butane", backend="HEOS", exceptions=True)


def _stagnation(fluid, p, q):
    st = fluid.get_state(jxp.PQ_INPUTS, p, q)
    return float(st.h), float(st.s)


def _synthetic(Th_in, Th_out, Tc_in, Tc_out, n, C_h, C_c):
    """Constant-cp counter-current profiles aligned cold inlet -> cold outlet."""
    Tc = np.linspace(Tc_in, Tc_out, n)
    Th = np.linspace(Th_out, Th_in, n)
    return ({"hot_side": {"states": {"T": Th, "h": Th * 1000.0}},
             "cold_side": {"states": {"T": Tc, "h": Tc * 1000.0}}}, C_h / 1000.0, C_c / 1000.0)


# ------------------------------------------------------------------ pump and turbine maps
def test_pump_map():
    p = cod.PumpMap(300.0, 0.01, 200.0, 1.2, 0.85)
    assert p.phi0 == pytest.approx(np.sqrt(6.0))
    assert p.efficiency(1.0) == pytest.approx(0.85, abs=1e-14)
    eps = 1e-6
    assert abs((p.efficiency(1 + eps) - p.efficiency(1 - eps)) / (2 * eps)) < 1e-6
    phi = np.linspace(1e-6, p.phi0 - 1e-6, 2001)
    assert np.all(p.efficiency(phi)[np.abs(phi - 1) > 1e-3] < 0.85)
    assert p.efficiency(0.0) == 0.0 and abs(p.efficiency(p.phi0)) < 1e-12
    r = np.array([0.5, 0.8, 1.0, 1.2])  # affinity at corresponding points
    H = p.head(0.9 * r * p.Qv_d, r * p.omega_d)
    assert np.allclose(H / r**2 / H[2], 1.0)
    assert cod.PumpMap.from_dict(p.to_dict()).to_dict() == p.to_dict()


def test_turbine_map():
    assert cod.turbine_efficiency(1.0, 0.9) == pytest.approx(0.9)
    assert cod.turbine_efficiency(0.0, 0.9) == 0.0 and cod.turbine_efficiency(2.0, 0.9) == 0.0
    assert cod.turbine_efficiency(0.9, 0.9) < 0.9 and cod.turbine_efficiency(1.1, 0.9) < 0.9


# ------------------------------------------------------------------ heat-exchanger conductance
def test_ua_matches_analytic_lmtd():
    prof, _, m_c = _synthetic(400.0, 340.0, 300.0, 360.0, 31, 2000.0, 2000.0)
    assert cod.integrate_ua(prof, m_c)["UA"] == pytest.approx(2000.0 * 60.0 / 40.0, rel=1e-12)
    prof, _, m_c = _synthetic(400.0, 320.0, 300.0, 340.0, 31, 1000.0, 2000.0)
    lmtd = ((400 - 340) - (320 - 300)) / np.log(60 / 20)
    assert cod.integrate_ua(prof, m_c)["UA"] == pytest.approx(1000.0 * 80.0 / lmtd, rel=1e-12)


def test_ua_continuous_and_monotone_through_crossing():
    ln_ua = []
    for th_in in np.linspace(361.0, 355.0, 61):  # hot-end approach 1 K -> -5 K
        prof, _, m_c = _synthetic(th_in, 320.0, 300.0, 360.0, 31, 1000.0, 1000.0)
        ln_ua.append(cod.integrate_ua(prof, m_c)["ln_UA"])
    steps = np.diff(ln_ua)
    assert np.all(np.isfinite(ln_ua)) and np.all(steps > 0)
    assert np.max(steps) < 1.0  # no jumps on a 0.1 K grid (slope about 1/delta per K)


def test_ua_finite_for_100_K_crossing():
    prof, _, m_c = _synthetic(260.0, 250.0, 300.0, 360.0, 31, 1000.0, 1000.0)  # approaches -40 to -100 K
    assert np.isfinite(cod.integrate_ua(prof, m_c)["ln_UA"])


# ------------------------------------------------------------------ nozzle
def _isentropic_flux(fluid, h0, s0, p):
    """Mass flux after isentropic expansion from (h0, s0) to p, computed independently."""
    st = fluid.get_state(jxp.PSmass_INPUTS, p, s0)
    return float(st.d) * np.sqrt(2.0 * (h0 - float(st.h)))


@pytest.mark.parametrize("name,p,q", NOZZLE_CASES)
def test_choked_flux_is_maximum_of_isentropic_flux(name, p, q):
    fl = jxp.Fluid(name=name, backend="HEOS", exceptions=True)
    h0, s0 = _stagnation(fl, p, q)
    flow = cod.nozzle_flow(fl, h0, s0, 0.05 * p, 200.0)
    G_star, p_star = flow["G"], flow["critical_pressure"]
    assert flow["choked"] and 0.05 * p < p_star < p
    # The flux at the critical pressure is the returned flux
    assert _isentropic_flux(fl, h0, s0, p_star) == pytest.approx(G_star, rel=1e-6)
    # No pressure gives a larger flux (fine grid around p*, coarse grid over the whole range)
    grid = np.concatenate([p_star * np.linspace(0.9, 1.1, 41), np.linspace(0.06 * p, 0.99 * p, 60)])
    assert max(_isentropic_flux(fl, h0, s0, pi) for pi in grid) <= G_star * (1 + 1e-6)


def test_unchoked_nozzle_flow(butane):
    h0, s0 = _stagnation(butane, 15e5, 0.05)
    choked = cod.nozzle_flow(butane, h0, s0, 1e5, 200.0)
    p_star, G_star = choked["critical_pressure"], choked["G"]
    above = cod.nozzle_flow(butane, h0, s0, p_star * (1 + 1e-4), 200.0)
    assert not above["choked"] and above["G"] == pytest.approx(G_star, rel=1e-4)
    assert cod.nozzle_flow(butane, h0, s0, 15e5 * 0.999, 200.0)["G"] < 0.2 * G_star
    G = [cod.nozzle_flow(butane, h0, s0, p, 200.0)["G"] for p in np.linspace(1.001 * p_star, 0.99 * 15e5, 10)]
    assert np.all(np.diff(G) < 0)


# ------------------------------------------------------------------ end-to-end example (slow)
@pytest.mark.skipif(os.environ.get("THERMOPT_RUN_SLOW") != "1", reason="set THERMOPT_RUN_SLOW=1 (about 3 min)")
def test_example_end_to_end(tmp_path, monkeypatch):
    for name in ("run_offdesign.py", "case_novec649_design.yaml", "case_novec649_offdesign.yaml"):
        shutil.copy(EXAMPLE / name, tmp_path)
    monkeypatch.chdir(tmp_path)
    ex = runpy.run_path("run_offdesign.py", run_name="example")  # design and sizing

    # Consistency: design conditions with the design controls fixed reproduce the design point
    cfg = ex["offdesign_config"](1.0, 0.0, 0.0)
    cfg["problem_formulation"]["fixed_parameters"]["controls"] = ex["design_controls"]
    for name in ex["design_controls"]:
        cfg["problem_formulation"]["design_variables"].pop(name)
    cycle = ex["offdesign"]
    cycle.load_config(cfg)
    cycle.run_optimization()
    assert cycle.solver.success
    W_net = cycle.problem.cycle_data["energy_analysis"]["net_system_power"]
    assert W_net == pytest.approx(ex["design_point"]["energy_analysis"]["net_system_power"], rel=1e-6)

    # Optimal controls at a different operating point
    expander = ex["f"](0.9, 5.0, 4.0)
    assert cycle.solver.success
    assert 0.02 < expander["inlet_quality"] < 0.98
