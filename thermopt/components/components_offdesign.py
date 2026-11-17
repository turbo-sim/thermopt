"""Off-design component models for fixed hardware: pump map, turbine map,
equivalent nozzle and heat-exchanger UA integration. SI units."""

import math

import jaxprop as jxp
import numpy as np
import pysolver_view as psv
from scipy.optimize import minimize_scalar

__all__ = [
    "GRAVITY",
    "CROSSING_THRESHOLD",
    "PumpMap",
    "turbine_efficiency",
    "find_choked_state",
    "nozzle_flow",
    "integrate_ua",
]

GRAVITY = 9.80665  # m/s^2
CROSSING_THRESHOLD = 0.2  # K, approach below which the crossing extension applies


# =============================================================================
# Pump
# =============================================================================
class PumpMap:
    r"""Variable-speed pump anchored at a design point (analytical, not validated against pump data).

    Similarity variables:
        r = omega / omega_d,  phi = Q_v / (r Q_v,d),  k = H_shutoff / H_d

    Head (parabola: k H_d at phi = 0, H_d at phi = 1, zero at phi0 = sqrt(k / (k - 1))):
        H = H_d r^2 [k - (k - 1) phi^2]

    Shaft power (assumed quadratic in phi):
        P = rho g Q_v,d H_d r^3 D(phi) / eta_d,  D(phi) = (k - 1/2) + (3/2 - k) phi^2
        D(1) = 1 gives P = rho g Q_v,d H_d / eta_d at design; the phi^2 coefficient
        puts the efficiency peak at phi = 1.

    Efficiency from eta = rho g Q_v H / P:
        eta = eta_d phi [k - (k - 1) phi^2] / D(phi)
        eta = eta_d at phi = 1, eta = 0 at phi = 0 and phi = phi0.

    Head is isentropic: H = [h(p2, s1) - h1] / g.
    """

    def __init__(self, design_speed, design_volumetric_flow, design_head,
                 shutoff_head_ratio, design_efficiency):
        self.omega_d = float(design_speed)
        self.Qv_d = float(design_volumetric_flow)
        self.H_d = float(design_head)
        self.k = float(shutoff_head_ratio)
        self.eta_d = float(design_efficiency)
        self.phi0 = math.sqrt(self.k / (self.k - 1.0))  # zero-head flow coefficient

    def to_dict(self):
        return {
            "design_speed": self.omega_d,
            "design_volumetric_flow": self.Qv_d,
            "design_head": self.H_d,
            "shutoff_head_ratio": self.k,
            "design_efficiency": self.eta_d,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(d["design_speed"], d["design_volumetric_flow"], d["design_head"],
                   d["shutoff_head_ratio"], d["design_efficiency"])

    def flow_coefficient(self, Qv, omega):
        return Qv / (omega / self.omega_d * self.Qv_d)

    def head(self, Qv, omega):
        phi = self.flow_coefficient(Qv, omega)
        return self.H_d * (omega / self.omega_d) ** 2 * (self.k - (self.k - 1.0) * phi**2)

    def efficiency(self, phi):
        k = self.k
        return self.eta_d * phi * (k - (k - 1.0) * phi**2) / ((k - 0.5) + (1.5 - k) * phi**2)


# =============================================================================
# Turbine
# =============================================================================
def turbine_efficiency(x, design_efficiency):
    """Velocity-ratio map eta = eta_d x (2 - x), with x = nu / nu_d the velocity ratio
    relative to design: peak eta_d at x = 1, zero at x = 0 and x = 2."""
    return design_efficiency * x * (2.0 - x)


# =============================================================================
# Equivalent nozzle
# =============================================================================
def find_choked_state(fluid, h0, s0, T_min):
    """Maximum isentropic mass flux for stagnation state (h0, s0).

    Returns ``(state_star, G_star, v_star)``. Adapted from
    ``flashing_nozzle_design.functions.find_choked_state`` (github.com/turbo-sim/
    flashing_nozzle_design, commit f2c3ceb): search up to v_max without a fallback interval,
    DomainError below T_min.
    """

    # Objective: maximize rho * v  ->  minimize its negative
    def neg_mass_flux(v):
        st = fluid.get_state(jxp.HmassSmass_INPUTS, h0 - 0.5 * v**2, s0)
        return -st.d * v

    # Check that the minimum temperature constraint is compatible
    st_min = fluid.get_state(jxp.SmassT_INPUTS, s0, T_min)
    if st_min.h > h0:
        raise psv.DomainError(f"Turbine inlet temperature below minimum ({T_min} K)")

    # Maximum possible velocity from energy balance
    v_max = np.sqrt(2 * (h0 - st_min.h))

    # Bounded search up to the velocity at T_min
    sol = minimize_scalar(neg_mass_flux, bounds=(1.0, v_max), method="bounded")

    v_star = sol.x
    st_star = fluid.get_state(jxp.HmassSmass_INPUTS, h0 - 0.5 * v_star**2, s0)

    return st_star, st_star.d * v_star, v_star,


def nozzle_flow(fluid, h0, s0, p_back, T_min):
    """Mass flux through the equivalent nozzle from stagnation state (h0, s0).

    G = G* (choked) if p_back <= p*, else the isentropic flux at p_back.
    """
    state_star, G_star, _ = find_choked_state(fluid, h0, s0, T_min)
    p_star = float(state_star.p)
    if p_back <= p_star:
        return {"G": float(G_star), "choked": True, "critical_pressure": p_star}
    st = fluid.get_state(jxp.PSmass_INPUTS, p_back, s0)
    dh = h0 - float(st.h)
    if dh <= 0.0:
        raise psv.DomainError(f"No isentropic enthalpy drop to {p_back:.4g} Pa")
    return {"G": float(st.d) * math.sqrt(2.0 * dh), "choked": False, "critical_pressure": p_star}


# =============================================================================
# Heat-exchanger conductance
# =============================================================================
def integrate_ua(hx, m_cold, delta=CROSSING_THRESHOLD):
    """Required UA of a ``heat_exchanger`` profile (cold inlet to cold outlet).

    UA = sum_i m_c (h_c,i+1 - h_c,i) / dT_lm,i. Approaches below ``delta`` are
    replaced by delta exp((dT - delta) / delta) (same value and slope at delta),
    so crossed profiles give a large but finite UA (overflow beyond about 700 delta
    of crossing); exact when all approaches are >= delta.
    """
    dT = np.asarray(hx["hot_side"]["states"]["T"], dtype=float) - np.asarray(hx["cold_side"]["states"]["T"], dtype=float)
    dT_ext = np.where(dT >= delta, dT, delta * np.exp((np.minimum(dT, delta) - delta) / delta))
    dQ = m_cold * np.diff(np.asarray(hx["cold_side"]["states"]["h"], dtype=float))
    a, b = dT_ext[:-1], dT_ext[1:]
    ratio = a / b
    equal = np.abs(ratio - 1.0) < 1e-6
    dT_lm = np.where(equal, 0.5 * (a + b), (a - b) / np.log(np.where(equal, 2.0, ratio)))
    UA = float(np.sum(dQ / dT_lm))
    return {"UA": UA, "ln_UA": math.log(UA), "min_approach": float(dT.min())}
