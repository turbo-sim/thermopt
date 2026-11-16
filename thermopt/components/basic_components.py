import numpy as np
import jaxprop as props

from scipy.integrate import solve_ivp
from ..utilities.fluid_utils import process_state_options
from .turbomachinery_nondimensional import solve_compressor, solve_turbine


def _snap_saturation_nodes(hot_side, cold_side, fluid_hot, fluid_cold, counter_current):
    """Move the wet node next to each phase boundary, then match the other stream.

    A single linear interpolation of h - h_sat across the original crossing
    estimates its pressure; a PQ evaluation puts the node exactly on saturation.
    Its enthalpy determines the common heat-duty coordinate. The opposite stream
    follows its original linear h/p path at that coordinate. With pressure loss,
    the saturated stream's pressure location is approximate (the saturation curve
    is nonlinear); there is no iterative solve. At constant pressure it is exact.

    Only original two-phase/single-phase transitions are considered. Endpoints
    and array lengths stay fixed. If boundaries compete for one wet node, retain
    the smaller temperature gap. A grid with no interior wet nodes cannot snap.
    """
    sides = (hot_side, cold_side)
    fluids = (fluid_hot, fluid_cold)
    original = [side["states"] for side in sides]
    count = len(original[0].enthalpy)
    dh = [float(side["state_out"].enthalpy - side["state_in"].enthalpy) for side in sides]
    if count < 3 or any(change == 0 for change in dh):
        return
    uniform = np.linspace(0, 1, count)  # Cold-inlet to cold-outlet coordinate.

    def fraction(side_index, x):
        return 1 - x if side_index == 0 and counter_current else x

    candidates = {}
    for side_index, (side, fluid, states) in enumerate(zip(sides, fluids, original)):
        if not hasattr(fluid, "triple_point_liquid"):
            continue
        pressure_min = float(fluid.triple_point_liquid.pressure)
        pressure_max = float(fluid.critical_point.pressure)
        two_phase = np.asarray(states.is_two_phase, dtype=bool)
        for left in range(count - 1):
            right = left + 1
            if two_phase[left] == two_phase[right]:
                continue
            current, outside = (left, right) if two_phase[left] else (right, left)
            if current in (0, count - 1) or not 0 < float(states.quality_mass[current]) < 1:
                continue  # Preserve endpoints and nodes already on saturation.
            p_wet, p_dry = float(states.pressure[current]), float(states.pressure[outside])
            if not (pressure_min <= p_wet < pressure_max
                    and pressure_min <= p_dry < pressure_max):
                continue

            # Single-phase trimmed Q is an entropy proxy in jaxprop, and can
            # label superheated dry fluids as liquid. Classify using enthalpy.
            sat_dry = fluid.get_state(props.PQ_INPUTS, p_dry, 0)
            quality = 0 if float(states.enthalpy[outside]) < float(sat_dry.enthalpy) else 1
            if quality == 1:
                sat_dry = fluid.get_state(props.PQ_INPUTS, p_dry, 1)
            saturated = sat_dry
            if p_wet != p_dry:
                sat_wet = fluid.get_state(props.PQ_INPUTS, p_wet, quality)
                r_dry = float(states.enthalpy[outside] - sat_dry.enthalpy)
                r_wet = float(states.enthalpy[current] - sat_wet.enthalpy)
                if r_wet == r_dry:
                    continue
                weight = -r_dry / (r_wet - r_dry)
                if not 0 < weight < 1:
                    continue
                p = p_dry + weight * (p_wet - p_dry)
                saturated = fluid.get_state(props.PQ_INPUTS, p, quality)

            x = fraction(side_index, (float(saturated.enthalpy) - float(side["state_in"].enthalpy)) / dh[side_index])
            if not uniform[left] < x < uniform[right]:
                continue  # Do not cross an existing node or duplicate an endpoint.
            other_index = 1 - side_index
            other = sides[other_index]
            other_fraction = fraction(other_index, x)
            other_h = float(other["state_in"].enthalpy) + other_fraction * dh[other_index]
            other_p = float(other["state_in"].pressure) + other_fraction * float(
                other["state_out"].pressure - other["state_in"].pressure)
            matched = fluids[other_index].get_state(props.HmassP_INPUTS, other_h, other_p)
            pair = (saturated, matched) if side_index == 0 else (matched, saturated)
            gap = float(pair[0].temperature - pair[1].temperature)
            if current not in candidates or gap < candidates[current][0]:
                candidates[current] = (gap, x, pair)
    if candidates:
        duty = uniform.copy()
        updated = [[states.at_index(i) for i in range(count)] for states in original]
        for index, (_, x, pair) in sorted(candidates.items()):
            if duty[index - 1] < x < duty[index + 1]:
                duty[index] = x
                updated[0][index], updated[1][index] = pair
        for side, states in zip(sides, updated):
            side["states"] = props.FluidState.stack(states)


def heat_exchanger(
    fluid_hot,
    h_in_hot,
    h_out_hot,
    p_in_hot,
    p_out_hot,
    fluid_cold,
    h_in_cold,
    h_out_cold,
    p_in_cold,
    p_out_cold,
    num_steps=50,
    counter_current=True,
    include_saturation_nodes=False,
):
    """
    Simulate a counter-current or co-current heat exchanger using discretized enthalpy and pressure profiles.

    The heat exchange is discretized along both streams. Thermodynamic states are evaluated at linearly spaced
    enthalpy and pressure points. The resulting states are ordered in the direction of increasing temperature.
    For counter-current flow, the hot-side state array is flipped to enable element-wise temperature difference
    computation.

    Parameters
    ----------
    fluid_hot : jaxprop.Fluid
        Fluid object for the hot side.
    h_in_hot : float
        Inlet specific enthalpy [J/kg] of the hot fluid.
    h_out_hot : float
        Outlet specific enthalpy [J/kg] of the hot fluid.
    p_in_hot : float
        Inlet pressure [Pa] of the hot fluid.
    p_out_hot : float
        Outlet pressure [Pa] of the hot fluid.
    fluid_cold : jaxprop.Fluid
        Fluid object for the cold side.
    h_in_cold : float
        Inlet specific enthalpy [J/kg] of the cold fluid.
    h_out_cold : float
        Outlet specific enthalpy [J/kg] of the cold fluid.
    p_in_cold : float
        Inlet pressure [Pa] of the cold fluid.
    p_out_cold : float
        Outlet pressure [Pa] of the cold fluid.
    num_steps : int, optional
        Number of discretization steps (default is 50).
    counter_current : bool, optional
        If True, assumes counter-current flow and flips hot-side state arrays (default is True).

    include_saturation_nodes : bool, optional
        False by default; the uniform calculation is unchanged when disabled.
        After uniform evaluation, snap the first interior two-phase node next
        to each single-phase region to Q=0 or Q=1. Estimate saturation pressure
        by linear interpolation across that crossing, then update the opposite
        stream's enthalpy AND pressure at the same heat-duty fraction.
        No root finding or extra nodes. Endpoints remain unchanged.
        Saturation is exact; its pressure location on the baseline path is
        approximate for nonzero pressure loss. Crossings without an interior
        two-phase node cannot snap. Conflicting candidates retain the smaller
        temperature gap, subject to heat-duty ordering. Retain enough nodes
        to resolve the two-phase regions and check accuracy with refinement.

    Returns
    -------
    dict
        Dictionary containing:
        - 'type': str, component type
        - 'hot_side': dict, state and property arrays for the hot stream
        - 'cold_side': dict, state and property arrays for the cold stream
        - 'q_hot_side': float, specific enthalpy drop [J/kg] on the hot side
        - 'q_cold_side': float, specific enthalpy gain [J/kg] on the cold side
        - 'temperature_hot_side': ndarray, temperatures [K] for the hot side (ordered by increasing T)
        - 'temperature_cold_side': ndarray, temperatures [K] for the cold side
        - 'temperature_difference': ndarray, element-wise temperature difference [K]
        - 'mass_flow_ratio': float, ratio of cold-side to hot-side mass flow required for heat balance
    """

    if include_saturation_nodes:
        for fluid in (fluid_hot, fluid_cold):
            if "saturation" not in fluid.capabilities.phase_operations:
                raise props.UnsupportedOperationError(
                    fluid_name=fluid.name, provider=fluid.provider, model=fluid.model,
                    request="include_saturation_nodes",
                    reason="Saturation-node snapping requires a fluid with saturation support.",
                    supported=("include_saturation_nodes=False", "source-EOS saturation"),
                )

    # Evaluate properties on the hot side
    hot_side = heat_transfer_process(
        fluid=fluid_hot,
        h_1=h_in_hot,
        h_2=h_out_hot,
        p_1=p_in_hot,
        p_2=p_out_hot,
        num_steps=num_steps,
    )

    # Evaluate properties in the cold side
    cold_side = heat_transfer_process(
        fluid=fluid_cold,
        h_1=h_in_cold,
        h_2=h_out_cold,
        p_1=p_in_cold,
        p_2=p_out_cold,
        num_steps=num_steps,
    )


    # Sort values for temperature difference calculation
    if counter_current:
        hot_side["states"] = hot_side["states"].flipped()

    if include_saturation_nodes:
        _snap_saturation_nodes(hot_side, cold_side, fluid_hot, fluid_cold, counter_current)

    # Compute temperature difference
    dT = hot_side["states"]["temperature"] - cold_side["states"]["temperature"]

    # Compute mass flow ratio
    # Use 1.0 for cases with no heat exchange (i.e., no recuperator)
    dh_hot = hot_side["state_in"].enthalpy - hot_side["state_out"].enthalpy
    dh_cold = cold_side["state_out"].enthalpy - cold_side["state_in"].enthalpy
    mass_ratio = 1.00 if (dh_cold == 0 or dh_hot == 0) else dh_cold / dh_hot

    # Create result dictionary
    result = {
        "type": "heat_exchanger",
        "hot_side": hot_side,
        "cold_side": cold_side,
        "q_hot_side": dh_hot,
        "q_cold_side": dh_cold,
        "temperature_hot_side": hot_side["states"]["temperature"],
        "temperature_cold_side": cold_side["states"]["temperature"],
        "temperature_difference": dT,
        "mass_flow_ratio": mass_ratio,
    }

    return result


def heat_transfer_process(fluid, h_1, p_1, h_2, p_2, num_steps=25):
    """
    Compute a discretized heat transfer process by idiscretizing enthalpy and pressure between inlet and outlet.

    This function generates a sequence of thermodynamic states along a heat transfer path
    by linearly spacing enthalpy and pressure between inlet and outlet values.
    The resulting states are organized in the direction of increasing enthalpy,
    which for sensible heating/cooling corresponds to increasing temperature.

    Parameters
    ----------
    fluid : jaxprop.Fluid
        jaxprop fluid object used for property evaluation.
    h_1 : float
        Inlet specific enthalpy [J/kg].
    p_1 : float
        Inlet pressure [Pa].
    h_2 : float
        Outlet specific enthalpy [J/kg].
    p_2 : float
        Outlet pressure [Pa].
    num_steps : int, optional
        Number of discretization steps (default is 25).

    Returns
    -------
    dict
        Dictionary containing:
        - 'states': dict of arrays with thermodynamic properties along the process
        - 'fluid_name': str, name of the fluid
        - 'state_in': jaxprop state object at the inlet
        - 'state_out': jaxprop state object at the outlet
        - 'mass_flow': float, set to NaN (to be populated externally)
        - 'heat_flow': float, set to NaN (to be populated externally)
        - 'color': str, set to 'black' (optional use for plotting)
    """

    # Generate linearly spaced arrays for pressure and enthalpy
    p_array = np.linspace(p_1, p_2, num_steps)
    h_array = np.linspace(h_1, h_2, num_steps)

    # Compute states for hot or cold side
    states = fluid.get_state(props.HmassP_INPUTS, h_array, p_array)

    # Create result dictionary
    result = {
        "states": states,
        "fluid_name": fluid.name,
        "state_in": states.at_index(0),
        "state_out": states.at_index(-1),
        "mass_flow": np.nan,
        "heat_flow": np.nan,
        "color": "black",
    }

    return result



def _validate_polytropic_steps(num_steps):
    """Require enough output nodes for both endpoints, once per component."""
    if isinstance(num_steps, (bool, np.bool_)) or not isinstance(num_steps, (int, np.integer)) or num_steps < 2:
        raise ValueError("Polytropic num_steps must be an integer greater than or equal to 2.")


def compression_process(
    fluid, h_in, p_in, p_out, efficiency=None, efficiency_type="isentropic", mass_flow=None, data_in=None, num_steps=10
):
    """Evaluate an adiabatic compression process through the supplied jaxprop fluid.

    h_in is specific enthalpy [J/kg]; p_in and p_out are total pressures [Pa].
    efficiency_type selects a constant isentropic efficiency, a constant
    polytropic efficiency integrated along pressure, or the non-dimensional
    radial-compressor correlation of Parisi et al. (2025), Section 3.
    For the first two models, efficiency is a fraction. The non-dimensional
    model instead requires mass_flow [kg/s] and accepts partial design
    parameters in data_in; it solves for outlet enthalpy as in Section 4.
    num_steps controls the polytropic path resolution (default 10).
    Polytropic paths require an integer num_steps >= 2. Equal inlet/outlet
    pressures return a constant full-state path without integration.

    Returns a component dictionary with inlet/outlet FluidStates, the state
    path, positive specific_work and isentropic_work [J/kg], pressure ratio,
    and design/performance data. The non-dimensional model's data_out also
    contains positive shaft power [W], efficiencies, size, speed and stress.
    The result's mass_flow field is assigned later by cycle energy accounting.
    Invalid efficiency selectors raise ValueError; the radial solver also
    validates flow, pressure direction and design inputs.
    """

    if efficiency_type == "polytropic":
        _validate_polytropic_steps(num_steps)
        # SciPy uses NumPy scalars. Keep one representation from the inlet
        # through its trial calls, avoiding extra weak-type JIT specializations.
        # Preserve dtype so JaxProp still rejects bool/string/complex inputs.
        h_in, p_in, p_out = (np.asarray(value) for value in (h_in, p_in, p_out))
    data_in = {} if data_in is None else data_in
    # Resolve optional reporting capabilities once, outside property calls.
    state_options = process_state_options(fluid, generalize_quality=(efficiency_type == "non-dimensional"))
    # Compute inlet state
    state_in = fluid.get_state(
        props.HmassP_INPUTS, h_in, p_in, **state_options,
    )
    if efficiency_type != "non-dimensional":
        # Only its enthalpy contributes to work; it is not a reporting state.
        state_out_is = fluid.get_state(props.PSmass_INPUTS, p_out, state_in.entropy)

    # Evaluate compression process
    if efficiency_type == "isentropic":
        h_out = state_in.enthalpy + (state_out_is.enthalpy - state_in.enthalpy) / efficiency
        state_out = fluid.get_state(props.HmassP_INPUTS, h_out, p_out, **state_options)
        states = state_in + state_out
        data_out = {"isentropic_efficiency": efficiency}

    elif efficiency_type == "polytropic":
        # Differential equation defining the polytropic compression
        def odefun(p, h):
            # The ODE uses density only. Request a complete core state without
            # optional saturation diagnostics on these discarded trial points.
            # SciPy's initial p is a Python float, later values are NumPy
            # scalars. Keep that distinction out of the property JIT cache.
            state = fluid.get_state(props.HmassP_INPUTS, h[0], np.float64(p))
            # The host solver needs a scalar derivative. Avoid launching tiny
            # device arithmetic operations between already synchronized calls.
            density = float(state.density)
            dhdp = 1.0 / (efficiency * density)
            return dhdp

        if p_in == p_out:
            pressures = np.full(num_steps, p_in)
            enthalpies = np.full(num_steps, state_in.enthalpy)
        else:
            # Solve polytropic compression differential equation
            sol = solve_ivp(
                odefun,
                [state_in.pressure, p_out],
                [state_in.enthalpy],
                t_eval=np.linspace(state_in.pressure, p_out, num_steps),
                method="RK45",
            )
            if not sol.success:
                raise RuntimeError(f"Polytropic compression integration failed: {sol.message}")
            pressures, enthalpies = sol.t, sol.y[0]

        # Reconstruct the complete reporting path in one native batch call.
        # Reporting diagnostics and quality conventions match the scalar path;
        # integration equations, RK45 tolerances and sampling are unchanged.
        states = fluid.get_state(props.HmassP_INPUTS, enthalpies, pressures, **state_options)
        state_in, state_out = states.at_index(0), states.at_index(-1)
        data_out = {"polytropic_efficiency": efficiency}
        
    elif efficiency_type == "non-dimensional":
        solution = solve_compressor(fluid, state_in, p_out, mass_flow, data_in)
        # Preserve cycle quality annotations without converting h -> T -> h.
        state_out = fluid.get_state(
            props.HmassP_INPUTS, solution["state_out"].enthalpy, p_out,
            **state_options,
        )
        state_out_is = solution["state_out_isentropic"]
        states = state_in + state_out
        data_out = solution["performance"]
        data_in = solution["design"]

    else:
        raise ValueError("Invalid efficiency_type. Use 'isentropic', 'polytropic' or 'non-dimensional'.")

    # Compute work
    isentropic_work = state_out_is.enthalpy - state_in.enthalpy
    specific_work = state_out.enthalpy - state_in.enthalpy

    # Create result dictionary
    result = {
        "type": "compressor",
        "fluid_name": fluid.name,
        "states": states,
        "state_in": state_in,
        "state_out": state_out,
        "efficiency": efficiency,
        "efficiency_type": efficiency_type,
        "specific_work": specific_work,
        "isentropic_work": isentropic_work,
        "pressure_ratio": state_out.pressure / state_in.pressure,
        "mass_flow": np.nan,
        "color": "black",
        "data_in": data_in,
        "data_out": data_in | data_out
    }

    return result


def expansion_process(
    fluid, h_in, p_in, p_out, efficiency=None, efficiency_type="isentropic", mass_flow=None, data_in=None, num_steps=50
):
    """Evaluate an adiabatic expansion process through the supplied jaxprop fluid.

    h_in is specific enthalpy [J/kg]; p_in and p_out are total pressures [Pa].
    efficiency_type selects a constant isentropic efficiency, a constant
    polytropic efficiency integrated along pressure, or the non-dimensional
    radial-turbine correlation of Parisi et al. (2025), Section 3.
    For the first two models, efficiency is a fraction. The non-dimensional
    model instead requires mass_flow [kg/s] and accepts partial design
    parameters in data_in; it solves for outlet enthalpy as in Section 4.
    num_steps controls the polytropic path resolution (default 50).
    Polytropic paths require an integer num_steps >= 2. Equal inlet/outlet
    pressures return a constant full-state path without integration.

    Returns a component dictionary with inlet/outlet FluidStates, the state
    path, positive specific_work and isentropic_work [J/kg], pressure ratio,
    and design/performance data. The non-dimensional model's data_out also
    contains positive shaft power [W], efficiencies, size, speed and stress.
    The result's mass_flow field is assigned later by cycle energy accounting.
    Invalid efficiency selectors raise ValueError; the radial solver also
    validates flow, pressure direction and design inputs.
    """
    if efficiency_type == "polytropic":
        _validate_polytropic_steps(num_steps)
        # Normalize once at the host-solver boundary without changing input
        # dtype validation or requesting a projected property API.
        h_in, p_in, p_out = (np.asarray(value) for value in (h_in, p_in, p_out))
    data_in = {} if data_in is None else data_in
    # Resolve optional reporting capabilities once, outside property calls.
    state_options = process_state_options(fluid, generalize_quality=True)
    # Compute inlet state
    state_in = fluid.get_state(props.HmassP_INPUTS, h_in, p_in, **state_options)
    if efficiency_type != "non-dimensional":
        # Only its enthalpy contributes to work; it is not a reporting state.
        state_out_is = fluid.get_state(props.PSmass_INPUTS, p_out, state_in.entropy)
    if efficiency_type == "isentropic":
        # Compute outlet state according to the definition of isentropic efficiency
        h_out = state_in.enthalpy - efficiency * (state_in.enthalpy - state_out_is.enthalpy)
        state_out = fluid.get_state(props.HmassP_INPUTS, h_out, p_out, **state_options)
        states = state_in + state_out
        data_out = {"isentropic_efficiency": efficiency}
        
    elif efficiency_type == "polytropic":
        # Differential equation defining the polytropic expansion
        def odefun(p, h):
            # The ODE uses density only; optional reporting diagnostics are
            # evaluated once on the returned path after successful integration.
            state = fluid.get_state(props.HmassP_INPUTS, h[0], np.float64(p))
            density = float(state.density)
            dhdp = efficiency / density
            return dhdp

        if p_in == p_out:
            pressures = np.full(num_steps, p_in)
            enthalpies = np.full(num_steps, state_in.enthalpy)
        else:
            # Solve polytropic expansion differential equation
            sol = solve_ivp(
                odefun,
                [state_in.pressure, p_out],
                [state_in.enthalpy],
                t_eval=np.linspace(state_in.pressure, p_out, num_steps),
                method="RK45",
            )
            if not sol.success:
                raise RuntimeError(f"Polytropic expansion integration failed: {sol.message}")
            pressures, enthalpies = sol.t, sol.y[0]

        # Preserve all reporting fields using the public full-state batch API.
        states = fluid.get_state(props.HmassP_INPUTS, enthalpies, pressures, **state_options)
        state_in, state_out = states.at_index(0), states.at_index(-1)
        data_out = {"polytropic_efficiency": efficiency}
        
    elif efficiency_type == "non-dimensional":
        solution = solve_turbine(fluid, state_in, p_out, mass_flow, data_in)
        # Preserve cycle quality annotations without converting h -> T -> h.
        state_out = fluid.get_state(
            props.HmassP_INPUTS, solution["state_out"].enthalpy, p_out,
            **state_options,
        )
        state_out_is = solution["state_out_isentropic"]
        states = state_in + state_out
        data_out = solution["performance"]
        data_in = solution["design"]
        
    else:
        raise ValueError("Invalid efficiency_type. Use 'isentropic', 'polytropic' or 'non-dimensional'.")

    # Compute work
    isentropic_work = state_in.enthalpy - state_out_is.enthalpy
    specific_work = state_in.enthalpy - state_out.enthalpy

    # Create result dictionary
    result = {
        "type": "expander",
        "fluid_name": fluid.name,
        "states": states,
        "state_in": state_in,
        "state_out": state_out,
        "efficiency": efficiency,
        "efficiency_type": efficiency_type,
        "specific_work": specific_work,
        "isentropic_work": isentropic_work,
        "pressure_ratio": state_in.pressure / state_out.pressure,
        "mass_flow": np.nan,
        "color": "black",
        "data_in": data_in,
        "data_out": data_in | data_out
    }

    return result


def isenthalpic_valve(state_in, p_out, fluid, N=50):
    """
    Simulate an isenthalpic throttling process across a valve.

    The enthalpy remains constant across the valve, and pressure is reduced from inlet to outlet.
    The function returns a list of thermodynamic states along the expansion

    This function will yield similar results as an expansion process with a polytropic/isentropic efficiency of 0 %.

    Parameters
    ----------
    state_in : jaxprop.FluidState
        Inlet thermodynamic state.
    p_out : float
        Outlet pressure [Pa].
    fluid : jaxprop.Fluid
        Fluid object for property evaluation.
    N : int, optional
        Number of pressure steps between inlet and outlet (default is 50).

    Returns
    -------
    jaxprop.FluidState
        List of thermodynamic states along the isenthalpic path.
    """  
    p_array = np.linspace(state_in.pressure, p_out, N)
    h_array = state_in.enthalpy + 0. * p_array  # TODO: should be fixed
    states = fluid.get_state(props.HmassP_INPUTS, h_array, p_array)
    return states


def postprocess_ode(t, y, ode_handle):
    """
    Reconstruct thermodynamic states from ODE solution vectors.

    This function applies a user-defined ODE handle to each time step to extract 
    the corresponding thermodynamic state based on the integrated solution. 
    It is useful for evaluating additional fluid properties not stored in the raw `y` vector.

    Parameters
    ----------
    t : ndarray
        Time array [s].
    y : ndarray
        Solution array of shape (n_states, n_time_steps).
    ode_handle : callable
        Function that takes (t_i, y_i) and returns a tuple (_, state), 
        where `state` is a jaxprop.FluidState.

    Returns
    -------
    jaxprop.FluidState
        Batched immutable state reconstructed from scalar ODE trajectory states.
    """ 
    # Collect individual states
    states = []
    for t_i, y_i in zip(t, y.T):
        _, state = ode_handle(t_i, y_i)
        states.append(state)

    # Combine into one batched FluidState
    return props.FluidState.stack(states)


def compute_component_energy_flows(components):
    """
    Compute power and heat flows for a set of components in a thermodynamic cycle.

    For heat exchangers, the function calculates local and total heat transfer on both hot and cold sides, 
    assuming precomputed enthalpy profiles. For turbomachinery (compressors or expanders), it evaluates 
    the mechanical power based on mass flow and specific work. The results are stored in-place in the 
    `components` dictionary.

    Parameters
    ----------
    components : dict
        Dictionary of components, each with a 'type' field and required properties such as mass flow, 
        enthalpy states, and specific work.

    Returns
    -------
    None
        The function modifies the `components` dictionary in place, adding:
        - 'power' [W]
        - 'heat_flow' and 'heat_flow_bis' [W]
        - 'heat_balance' [W] for heat exchangers
    """
    for name, component in components.items():
        if component["type"] == "heat_exchanger":
            hot = component["hot_side"]
            cold = component["cold_side"]
            Q_hot_array = hot["mass_flow"] * (hot["state_in"]["enthalpy"] - hot["states"]["enthalpy"])
            Q_cold_array = cold["mass_flow"] * (cold["states"]["enthalpy"] - cold["state_in"]["enthalpy"])
            component["hot_side"]["heat_flow"] = Q_hot_array
            component["cold_side"]["heat_flow"] = Q_cold_array
            # component["hot_side"]["states"]["heat_flow"] = Q_hot_array
            # component["cold_side"]["states"]["heat_flow"] = Q_cold_array
            component["heat_flow"] = Q_hot_array[0]
            component["heat_flow_bis"] = Q_cold_array[-1]
            component["heat_balance"] = Q_hot_array[0] - Q_cold_array[-1]
            component["power"] = 0.0

        elif component["type"] in ["compressor", "expander"]:
            component["power"] = component["mass_flow"] * component["specific_work"]
            component["heat_flow_rate"] = 0.0
