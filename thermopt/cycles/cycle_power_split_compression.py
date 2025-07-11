import copy
import coolpropx as cpx

from .. import utilities

from ..components import compression_process, expansion_process, heat_exchanger, compute_component_energy_flows

COLORS_MATLAB = utilities.COLORS_MATLAB


def evaluate_cycle(
    variables,
    parameters,
    constraints,
    objective_function,
):
    # Create copies to not change the originals
    variables = copy.deepcopy(variables)
    parameters = copy.deepcopy(parameters)

    # Initialize fluid objects
    working_fluid = cpx.Fluid(
        **parameters.pop("working_fluid"), identifier="working_fluid"
    )
    heating_fluid = cpx.Fluid(
        **parameters.pop("heating_fluid"), identifier="heating_fluid"
    )
    cooling_fluid = cpx.Fluid(
        **parameters.pop("cooling_fluid"), identifier="cooling_fluid"
    )

    # Extract special points
    special_points = parameters.pop("special_points")

    # Extract heat source/sink parameters and give short names
    T_source_out_min = parameters["heat_source"].pop("minimum_temperature")
    p_source_out = parameters["heat_source"].pop("exit_pressure")
    p_sink_out = parameters["heat_sink"].pop("exit_pressure")

    # Compute coldest state at the heat source exit
    source_out_min = heating_fluid.get_state(
        cpx.PT_INPUTS, p_source_out, T_source_out_min
    )

    # Extract pressure drops and give short names
    dp_heater_h = parameters["heater"].pop("pressure_drop_hot_side")
    dp_heater_c = parameters["heater"].pop("pressure_drop_cold_side")
    dp_recup_lowT_h = parameters["recuperator_lowT"].pop("pressure_drop_hot_side")
    dp_recup_lowT_c = parameters["recuperator_lowT"].pop("pressure_drop_cold_side")
    dp_recup_highT_h = parameters["recuperator_highT"].pop("pressure_drop_hot_side")
    dp_recup_highT_c = parameters["recuperator_highT"].pop("pressure_drop_cold_side")
    dp_cooler_h = parameters["cooler"].pop("pressure_drop_hot_side")
    dp_cooler_c = parameters["cooler"].pop("pressure_drop_cold_side")

    # Extract design variables from dictionary (make sure all are used)
    expander_inlet_p = variables.pop("expander_inlet_pressure")
    expander_inlet_h = variables.pop("expander_inlet_enthalpy")
    main_compressor_inlet_p = variables.pop("main_compressor_inlet_pressure")
    main_compressor_inlet_h = variables.pop("main_compressor_inlet_enthalpy")
    split_compressor_inlet_h = variables.pop("split_compressor_inlet_enthalpy")
    recup_intermediate_h = variables.pop("recuperator_intermediate_enthalpy")
    mass_split_fraction = variables.pop("mass_split_fraction")
    heat_source_temperature_out = variables.pop("heat_source_exit_temperature")
    heat_sink_temperature_out = variables.pop("heat_sink_exit_temperature")

    # Evaluate main compressor
    dp = (1.0 - dp_heater_c) * (1.0 - dp_recup_highT_c) * (1.0 - dp_recup_lowT_c)
    main_compressor_outlet_p = expander_inlet_p / dp
    main_compressor_eff = parameters["main_compressor"].pop("efficiency")
    main_compressor_eff_type = parameters["main_compressor"].pop("efficiency_type")
    main_compressor = compression_process(
        working_fluid,
        main_compressor_inlet_h,
        main_compressor_inlet_p,
        main_compressor_outlet_p,
        main_compressor_eff,
        efficiency_type=main_compressor_eff_type,
    )

    # Evaluate re-compressor
    split_compressor_inlet_p = main_compressor_inlet_p / (1 - dp_cooler_h)
    dp = (1.0 - dp_heater_c) * (1.0 - dp_recup_highT_c)
    split_compressor_outlet_p = expander_inlet_p / dp
    split_compressor_eff = parameters["split_compressor"].pop("efficiency")
    split_compressor_eff_type = parameters["split_compressor"].pop("efficiency_type")
    split_compressor = compression_process(
        working_fluid,
        split_compressor_inlet_h,
        split_compressor_inlet_p,
        split_compressor_outlet_p,
        split_compressor_eff,
        split_compressor_eff_type,
    )

    # Evaluate expander
    dp = (1.0 - dp_cooler_h) * (1.0 - dp_recup_lowT_h) * (1.0 - dp_recup_highT_h)
    expander_outlet_p = main_compressor_inlet_p / dp
    expander_efficiency = parameters["expander"].pop("efficiency")
    expander_efficiency_type = parameters["expander"].pop("efficiency_type")
    expander = expansion_process(
        working_fluid,
        expander_inlet_h,
        expander_inlet_p,
        expander_outlet_p,
        expander_efficiency,
        expander_efficiency_type,
    )

    # Evaluate low temperature recuperator
    h_in_hot = recup_intermediate_h
    p_in_hot = expander_outlet_p * (1 - dp_recup_highT_h)
    h_out_hot = split_compressor_inlet_h
    p_out_hot = split_compressor_inlet_p
    h_in_cold = main_compressor["state_out"].h
    p_in_cold = main_compressor["state_out"].p
    h_out_cold = h_in_cold + (1 / (1 - mass_split_fraction)) * (h_in_hot - h_out_hot)
    p_out_cold = split_compressor_outlet_p
    num_elements = parameters["recuperator_lowT"].pop("num_elements")
    recuperator_lowT = heat_exchanger(
        working_fluid,
        h_in_hot,
        h_out_hot,
        p_in_hot,
        p_out_hot,
        working_fluid,
        h_in_cold,
        h_out_cold,
        p_in_cold,
        p_out_cold,
        counter_current=True,
        num_steps=num_elements,
    )

    # Evaluate high temperature recuperator
    h_in_cold = (
        mass_split_fraction * split_compressor["state_out"].h
        + (1 - mass_split_fraction) * recuperator_lowT["cold_side"]["state_out"].h
    )
    h_out_cold = (expander["state_out"].h - recup_intermediate_h) + h_in_cold
    p_in_cold = split_compressor["state_out"].p
    p_out_cold = expander_inlet_p / (1 - dp_heater_c)
    h_in_hot = expander["state_out"].h
    h_out_hot = recup_intermediate_h
    p_in_hot = expander_outlet_p
    p_out_hot = p_in_hot * (1 - dp_recup_highT_h)
    num_elements = parameters["recuperator_highT"].pop("num_elements")
    recuperator_highT = heat_exchanger(
        working_fluid,
        h_in_hot,
        h_out_hot,
        p_in_hot,
        p_out_hot,
        working_fluid,
        h_in_cold,
        h_out_cold,
        p_in_cold,
        p_out_cold,
        counter_current=True,
        num_steps=num_elements,
    )

    # Evaluate heater
    h_in_cold = recuperator_highT["cold_side"]["state_out"].h
    p_in_cold = recuperator_highT["cold_side"]["state_out"].p
    h_out_cold = expander["state_in"].h
    p_out_cold = expander["state_in"].p
    T_in_hot = parameters["heat_source"].pop("inlet_temperature")
    p_in_hot = parameters["heat_source"].pop("inlet_pressure")
    h_in_hot = heating_fluid.get_state(cpx.PT_INPUTS, p_in_hot, T_in_hot).h
    T_out_hot = heat_source_temperature_out
    p_out_hot = p_in_hot * (1 - dp_heater_h)
    h_out_hot = heating_fluid.get_state(cpx.PT_INPUTS, p_out_hot, T_out_hot).h
    num_elements = parameters["heater"].pop("num_elements")
    heater = heat_exchanger(
        heating_fluid,
        h_in_hot,
        h_out_hot,
        p_in_hot,
        p_out_hot,
        working_fluid,
        h_in_cold,
        h_out_cold,
        p_in_cold,
        p_out_cold,
        counter_current=True,
        num_steps=num_elements,
    )

    # Evaluate heat source pump
    h_in = heater["hot_side"]["state_out"].h
    p_in = heater["hot_side"]["state_out"].p
    p_out = p_source_out
    efficiency = parameters["heat_source_pump"].pop("efficiency")
    efficiency_type = parameters["heat_source_pump"].pop("efficiency_type")
    heat_source_pump = compression_process(
        heating_fluid,
        h_in,
        p_in,
        p_out,
        efficiency,
        efficiency_type,
    )

    # Evaluate heat sink pump
    T_in = parameters["heat_sink"].pop("inlet_temperature")
    p_in = parameters["heat_sink"].pop("inlet_pressure")
    h_in = cooling_fluid.get_state(cpx.PT_INPUTS, p_in, T_in).h
    p_out = p_sink_out / (1 - dp_cooler_c)
    efficiency = parameters["heat_sink_pump"].pop("efficiency")
    efficiency_type = parameters["heat_sink_pump"].pop("efficiency_type")
    heat_sink_pump = compression_process(
        cooling_fluid,
        h_in,
        p_in,
        p_out,
        efficiency,
        efficiency_type,
    )

    # Evaluate cooler
    p_in_cold = heat_sink_pump["state_out"].p
    h_in_cold = heat_sink_pump["state_out"].h
    p_out_cold = p_in_cold * (1 - dp_cooler_c)
    T_out_cold = heat_sink_temperature_out
    h_out_cold = cooling_fluid.get_state(cpx.PT_INPUTS, p_out_cold, T_out_cold).h
    h_in_hot = recuperator_lowT["hot_side"]["state_out"].h
    p_in_hot = recuperator_lowT["hot_side"]["state_out"].p
    h_out_hot = main_compressor["state_in"].h
    p_out_hot = main_compressor["state_in"].p
    num_elements = parameters["cooler"].pop("num_elements")
    cooler = heat_exchanger(
        working_fluid,
        h_in_hot,
        h_out_hot,
        p_in_hot,
        p_out_hot,
        cooling_fluid,
        h_in_cold,
        h_out_cold,
        p_in_cold,
        p_out_cold,
        counter_current=True,
        num_steps=num_elements,
    )

    # Compute mass flow rates
    W_net = parameters.pop("net_power")
    expander_work = expander["specific_work"]
    main_compression_work = main_compressor["specific_work"] * (1 - mass_split_fraction)
    split_compression_work = split_compressor["specific_work"] * (mass_split_fraction)
    m_total = W_net / (expander_work - main_compression_work - split_compression_work)
    m_main = m_total * (1 - mass_split_fraction)
    m_split = m_total * mass_split_fraction
    m_source = m_total * heater["mass_flow_ratio"]
    m_sink = m_main / cooler["mass_flow_ratio"]

    # Add the mass flow to the components
    heater["hot_side"]["mass_flow"] = m_source
    heater["cold_side"]["mass_flow"] = m_total
    recuperator_lowT["hot_side"]["mass_flow"] = m_total
    recuperator_lowT["cold_side"]["mass_flow"] = m_main
    recuperator_highT["hot_side"]["mass_flow"] = m_total
    recuperator_highT["cold_side"]["mass_flow"] = m_total
    cooler["hot_side"]["mass_flow"] = m_main
    cooler["cold_side"]["mass_flow"] = m_sink
    expander["mass_flow"] = m_total
    main_compressor["mass_flow"] = m_main
    split_compressor["mass_flow"] = m_split
    heat_source_pump["mass_flow"] = m_source
    heat_sink_pump["mass_flow"] = m_sink

    # Summary of components
    components = {
        "expander": expander,
        "split_compressor": split_compressor,
        "main_compressor": main_compressor,
        "recuperator_lowT": recuperator_lowT,
        "recuperator_highT": recuperator_highT,
        "heater": heater,
        "cooler": cooler,
        "heat_source_pump": heat_source_pump,
        "heat_sink_pump": heat_sink_pump,
    }
    compute_component_energy_flows(components)

    # First-law analysis
    Q_in = heater["heat_flow"]
    Q_out = cooler["heat_flow"]
    W_out = expander["power"]
    W_comp = main_compressor["power"] + split_compressor["power"]
    W_aux = heat_source_pump["power"] + heat_sink_pump["power"]
    W_in = W_comp + W_aux
    Q_in_max = m_source * (heater["hot_side"]["state_in"].h - source_out_min.h)
    cycle_efficiency = (W_out - W_in) / Q_in
    system_efficiency = (W_out - W_in) / Q_in_max
    backwork_ratio = W_comp / W_out
    energy_balance = (Q_in + W_in) - (W_out + Q_out)

    # Mixing chamber check
    # TODO remember to include mixing chamber exergy destruction in Second Law analysis
    a = m_total * recuperator_highT["cold_side"]["state_in"].h
    b = m_main * recuperator_lowT["cold_side"]["state_out"].h
    c = m_split * split_compressor["state_out"].h
    error = a - b - c
    if abs(error) > 1e-4:
        raise ValueError(
            f"The mixing chamber energy balance does not close (error: {error:0.3}). Check implementation."
        )

    # Define dictionary with 1st Law analysis
    energy_analysis = {
        "heater_heat_flow": heater["heat_flow"],
        "heater_heat_flow_max": Q_in_max,
        "cooler_heat_flow": cooler["heat_flow"],
        "expander_power": expander["power"],
        "main_compressor_power": main_compressor["power"],
        "split_compressor_power": split_compressor["power"],
        "heat_source_pump_power": heat_source_pump["power"],
        "heat_sink_pump_power": heat_sink_pump["power"],
        "net_cycle_power": W_net,
        "net_system_power": W_out - W_in,
        "mass_flow_heat_source": m_source,
        "mass_flow_heat_sink": m_sink,
        "mass_flow_total": m_total,
        "mass_flow_main": m_main,
        "mass_flow_split": m_split,
        "split_fraction": mass_split_fraction,
        "cycle_efficiency": cycle_efficiency,
        "system_efficiency": system_efficiency,
        "backwork_ratio": backwork_ratio,
        "energy_balance": energy_balance,
    }

    output = {
        "components": components,
        "energy_analysis": energy_analysis,
    }

    # Evaluate objective function and constraints
    output = {"components": components, "energy_analysis": energy_analysis}
    f = utilities.evaluate_objective_function(output, objective_function)
    c_eq, c_ineq, constraint_report = utilities.evaluate_constraints(output, constraints)


    # Set colors for plotting
    orange = COLORS_MATLAB[1]
    blue = COLORS_MATLAB[0]
    red = COLORS_MATLAB[6]
    heater["hot_side"]["plot_params"] = {"color": red, "linestyle": "-"}
    heater["cold_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    recuperator_lowT["hot_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    recuperator_lowT["cold_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    recuperator_highT["hot_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    recuperator_highT["cold_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    cooler["hot_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    cooler["cold_side"]["plot_params"] = {"color": blue, "linestyle": "-"}
    expander["plot_params"] = {"color": orange, "linestyle": "-"}
    main_compressor["plot_params"] = {"color": orange, "linestyle": "-"}
    split_compressor["plot_params"] = {"color": orange, "linestyle": "-"}

    # Check if any fixed parameter or design variable was not used
    utilities.check_for_unused_keys(parameters, "parameters", raise_error=True)
    utilities.check_for_unused_keys(variables, "variables", raise_error=True)

    # Cycle performance summary
    output = {
        **output,
        "working_fluid": working_fluid,
        "heating_fluid": heating_fluid,
        "cooling_fluid": cooling_fluid,
        "components": components,
        "objective_function": f,
        "equality_constraints": c_eq,
        "inequality_constraints": c_ineq,
        "constraints_report": constraint_report,
    }
    
    return output
