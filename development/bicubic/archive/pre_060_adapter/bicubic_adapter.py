"""Scoped JAXProp bicubic integration for the recuperated-cycle comparison."""

from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import jaxprop as jxp
import numpy as np

from thermopt.cycles import cycle_power_recuperated


DEFAULT_TABLE_DIR = Path(__file__).resolve().parent / "results" / "tables"
CORE_PROPERTIES = ("temperature", "entropy", "density")
CO2_DOMAIN = dict(h_min=110e3, h_max=1.3e6, p_min=8.0e6, p_max=39.0e6)
INPUT_PAIRS = {
    jxp.HmassP_INPUTS: ("enthalpy", "pressure"),
    jxp.PT_INPUTS: ("pressure", "temperature"),
    jxp.PSmass_INPUTS: ("pressure", "entropy"),
}


def _record_range(ranges, name, values):
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Nonfinite {name} in bicubic property query")
    limits = ranges.setdefault(name, {"minimum": float("inf"), "maximum": -float("inf")})
    limits["minimum"] = min(limits["minimum"], float(np.min(values)))
    limits["maximum"] = max(limits["maximum"], float(np.max(values)))


def _check_domain(model, enthalpy, pressure):
    for name, values, lower, upper in (
        ("enthalpy", enthalpy, model.h_min, model.h_max),
        ("pressure", pressure, model.p_min, model.p_max),
    ):
        values = np.asarray(values, dtype=float)
        tolerance = 1e-9 * max(abs(lower), abs(upper))
        if not np.all(np.isfinite(values)) or np.any(values < lower - tolerance) or np.any(values > upper + tolerance):
            raise ValueError(
                f"{model.fluid_name} bicubic {name} outside table domain "
                f"[{lower:.12g}, {upper:.12g}]. "
                f"Query range [{np.min(values):.12g}, {np.max(values):.12g}]"
            )


def _make_model(fluid_name, domain, grid_shape, table_dir):
    """Collect generation failures before accepting a generated or cached table."""
    generation_failures = []
    generation_count = [0]

    class AuditedFluidBicubic(jxp.FluidBicubic):
        @staticmethod
        def _gradients_central(fluid, h, p, eps_h, eps_p):
            try:
                gradients = jxp.FluidBicubic._gradients_central(fluid, h, p, eps_h, eps_p)
                for name in CORE_PROPERTIES:
                    if not np.all(np.isfinite(np.asarray(gradients[name], dtype=float))):
                        raise ValueError(f"Nonfinite {name} value or derivative")
            except Exception as error:
                generation_failures.append(dict(h=float(h), p=float(p), error=str(error)))
                raise
            generation_count[0] += 1
            return gradients

    n_h, n_p = grid_shape
    table_name = f"comparison_v1_{fluid_name}_HEOS_central_{n_h}x{n_p}"
    table_path = table_dir / f"{table_name}.pkl"
    loaded_from_cache = table_path.exists()
    model = AuditedFluidBicubic(
        fluid_name=fluid_name,
        backend="HEOS",
        **domain,
        N_h=n_h,
        N_p=n_p,
        gradient_method="central",
        table_name=table_name,
        table_dir=str(table_dir),
    )
    if generation_failures:
        table_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"Rejected {fluid_name} bicubic table with {len(generation_failures)} "
            f"failed nodes. First failure {generation_failures[0]}"
        )
    for name in CORE_PROPERTIES:
        for field in ("value", "grad_h", "grad_p", "grad_logP", "grad_ph", "grad_hlogP", "coeffs"):
            values = np.asarray(model.table[name][field])
            if not np.all(np.isfinite(values)):
                table_path.unlink(missing_ok=True)
                raise RuntimeError(f"Rejected {fluid_name} table with nonfinite {name}/{field}")
    return model, {
        "fluid_name": fluid_name,
        "implementation": "jaxprop.FluidBicubic",
        "table_path": str(table_path),
        "loaded_from_cache": loaded_from_cache,
        "grid_shape": [n_h, n_p],
        "domain": dict(domain),
        "gradient_method": "central",
        "generation_successful_nodes": generation_count[0],
        "generation_failed_nodes": len(generation_failures),
        "validated_properties": list(CORE_PROPERTIES),
        "query_count": 0,
        "state_count": 0,
        "query_counts_by_input_pair": {},
        "requested_ranges": {},
        "returned_ranges": {},
        "maximum_input_residuals": {},
        "query_failures": [],
    }


class BicubicFluidAdapter:
    """Preserve ThermOpt's fluid interface while auditing native interpolation."""

    def __init__(self, model, identifier, audit):
        self.model = model
        self.name = model.fluid_name
        self.backend = "HEOS"
        self.identifier = identifier if identifier is not None else self.name
        self.audit = audit

    def get_state(self, input_type, value_1, value_2, generalize_quality=False, supersaturation=False):
        if input_type not in INPUT_PAIRS:
            raise ValueError(f"Unsupported bicubic comparison input pair {input_type}")
        names = INPUT_PAIRS[input_type]
        pair_name = "-".join(names)
        requested = dict(zip(names, np.broadcast_arrays(value_1, value_2)))
        self.audit["query_count"] += 1
        self.audit["state_count"] += int(np.size(requested[names[0]]))
        pair_counts = self.audit["query_counts_by_input_pair"]
        pair_counts[pair_name] = pair_counts.get(pair_name, 0) + 1
        try:
            for name, values in requested.items():
                _record_range(self.audit["requested_ranges"], name, values)
            pressure = requested["pressure"]
            enthalpy = requested.get("enthalpy", self.model.h_min)
            _check_domain(self.model, enthalpy, pressure)
            state = self.model.get_state(input_type, value_1, value_2)
            for name in ("enthalpy", "pressure", *CORE_PROPERTIES):
                _record_range(self.audit["returned_ranges"], name, state[name])
            _check_domain(self.model, state.h, state.p)
            for name, values in requested.items():
                residual = np.abs(np.asarray(state[name]) - values)
                maximum_residual = float(np.max(residual))
                previous = self.audit["maximum_input_residuals"].get(name, 0.0)
                self.audit["maximum_input_residuals"][name] = max(previous, maximum_residual)
                tolerance = 1e-7 + 1e-9 * np.abs(values)
                if np.any(residual > tolerance):
                    raise ValueError(
                        f"{self.name} bicubic {pair_name} input recovery failed "
                        f"for {name}. Maximum absolute residual {maximum_residual:.12g}"
                    )
        except Exception as error:
            self.audit["query_failures"].append(dict(input_pair=pair_name, error=str(error)))
            raise
        return jxp.FluidState(
            fluid_name=self.name,
            identifier=self.identifier,
            **state.to_dict(),
        )


@contextmanager
def bicubic_cycle(grid_size=100, table_dir=DEFAULT_TABLE_DIR, all_fluids=False):
    """Use shared bicubic tables inside recuperated-cycle evaluations.

    ThermOpt's optimizer setup keeps HEOS special points, initial values and
    bounds. Native FluidBicubic PH, PS and PT evaluation supplies the selected
    cycle streams. Saturation and generalized-quality fields are unused by
    this example's objective and constraints.
    """
    grid_size = int(grid_size)
    if grid_size < 2:
        raise ValueError("Bicubic grid_size must be at least two")
    table_dir = Path(table_dir).resolve()
    table_dir.mkdir(parents=True, exist_ok=True)
    specifications = {"CO2": (CO2_DOMAIN, (grid_size, grid_size))}
    if all_fluids:
        specifications.update({
            "nitrogen": (dict(h_min=250e3, h_max=1e6, p_min=90e3, p_max=115e3), (60, 8)),
            "water": (dict(h_min=60e3, h_max=180e3, p_min=90e3, p_max=115e3), (60, 8)),
        })
    models = {}
    audit = {
        "implementation": "jaxprop.FluidBicubic",
        "all_fluids": bool(all_fluids),
        "setup_special_points_backend": "HEOS",
        "native_inverse_solvers": True,
        "extra_state_fields_requested_by_thermopt": ["generalize_quality", "supersaturation"],
        "extra_state_fields_evaluated": False,
        "fluids": {},
    }
    for name, (domain, shape) in specifications.items():
        model, model_audit = _make_model(name, domain, shape, table_dir)
        models[name.casefold()] = model
        audit["fluids"][name] = model_audit
    original_namespace = cycle_power_recuperated.cpx

    def fluid_factory(name, backend="HEOS", exceptions=True, identifier=None):
        model = models.get(name.casefold())
        if model is None:
            return original_namespace.Fluid(
                name=name, backend=backend, exceptions=exceptions, identifier=identifier
            )
        if backend != "HEOS":
            raise ValueError(f"Bicubic comparison requires HEOS tables, received {backend}")
        return BicubicFluidAdapter(model, identifier, audit["fluids"][model.fluid_name])

    cycle_power_recuperated.cpx = SimpleNamespace(
        Fluid=fluid_factory, PT_INPUTS=original_namespace.PT_INPUTS
    )
    try:
        yield audit
    finally:
        cycle_power_recuperated.cpx = original_namespace
