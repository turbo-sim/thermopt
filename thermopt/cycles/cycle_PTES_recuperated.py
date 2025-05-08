import copy
from .. import utilities
from .. import properties as props

from .components import compression_process, expansion_process, heat_exchanger

COLORS_MATLAB = utilities.COLORS_MATLAB


def evaluate_cycle(
    variables,
    parameters,
    constraints,
    objective_function,
    recuperated=True,
):
    # Create copies to not change the originals
    variables = copy.deepcopy(variables)
    parameters = copy.deepcopy(parameters)

    #print("hello")

    # Initialize fluid objects
    working_fluid = props.Fluid(
        **parameters.pop("working_fluid"), identifier="working_fluid"
    )
    heating_fluid = props.Fluid(
        **parameters.pop("heating_fluid"), identifier="heating_fluid"
    )
    cooling_fluid = props.Fluid(
        **parameters.pop("cooling_fluid"), identifier="cooling_fluid"
    )

    # Extract heat source/sink parameters and give short names
    # T_source_out_max = parameters["heat_source"].pop("exit_temperature_max")
    # T_source_out_min = parameters["heat_source"].pop("exit_temperature_min")
    # T_sink_out_max = parameters["heat_sink"].pop("exit_temperature_max")
    # T_sink_out_min = parameters["heat_sink"].pop("exit_temperature_min")
    p_source_out = parameters["heat_source"].pop("exit_pressure")
    p_sink_out = parameters["heat_sink"].pop("exit_pressure")

    # Compute coldest state at the heat source exit
    # source_out_min = heating_fluid.get_state(
    #     props.PT_INPUTS, p_source_out, T_source_out_min
    # )

    # Extract pressure drops and give short names
    dp_heater_h = parameters["heater_charge"].pop("pressure_drop_hot_side")
    dp_heater_c = parameters["heater_charge"].pop("pressure_drop_cold_side")
    dp_cooler_h = parameters["cooler_charge"].pop("pressure_drop_hot_side")
    dp_cooler_c = parameters["cooler_charge"].pop("pressure_drop_cold_side")
    if recuperated:
        dp_recup_h = parameters["recuperator_charge"].pop("pressure_drop_hot_side")
        dp_recup_c = parameters["recuperator_charge"].pop("pressure_drop_cold_side")
    else:
        dp_recup_c, dp_recup_h = 0.0, 0.0

    # Extract design variables from dictionary (make sure all are used)
    expander_inlet_p = variables.pop("expander_inlet_pressure_charge")
    expander_inlet_h = variables.pop("expander_inlet_enthalpy_charge")
    compressor_inlet_p = variables.pop("compressor_inlet_pressure_charge")
    compressor_inlet_h = variables.pop("compressor_inlet_enthalpy_charge")
    heat_source_temperature_out = variables.pop("heat_source_exit_temperature_charge")
    heat_sink_temperature_out = variables.pop("heat_sink_exit_temperature_charge")
    if recuperated:
        recuperator_effectiveness = variables.pop("recuperator_effectiveness_charge")
    else:
        recuperator_effectiveness = 0.0

    # Evaluate compressor
    dp = (1.0 - dp_cooler_h) * (1.0 - dp_recup_h)
    compressor_outlet_p = expander_inlet_p / dp
    compressor_eff = parameters["compressor_charge"].pop("efficiency")
    compressor_eff_type = parameters["compressor_charge"].pop("efficiency_type")
    compressor_charge = compression_process(
        working_fluid,
        compressor_inlet_h,
        compressor_inlet_p,
        compressor_outlet_p,
        compressor_eff,
        compressor_eff_type,
    )
    # utilities.print_dict(compressor["state_out"])

    # Evaluate expander
    dp = (1.0 - dp_heater_c) * (1.0 - dp_recup_c)
    expander_outlet_p = compressor_inlet_p / dp
    expander_efficiency = parameters["expander_charge"].pop("efficiency")
    expander_efficiency_type = parameters["expander_charge"].pop("efficiency_type")
    expander_charge = expansion_process(
        working_fluid,
        expander_inlet_h,
        expander_inlet_p,
        expander_outlet_p,
        expander_efficiency,
        expander_efficiency_type,
    )

    # Evaluate recuperator
    eps = recuperator_effectiveness
    T_out_hot = expander_charge["state_in"].T
    p_out_cold = compressor_charge["state_in"].p
    h_out_cold = compressor_charge["state_in"].h
    p_in_cold = p_out_cold / (1.0 - dp_recup_c)
    h_in_cold_ideal = working_fluid.get_state(props.PT_INPUTS, p_in_cold, T_out_hot).h
    h_in_cold_actual = h_out_cold - eps * (h_out_cold - h_in_cold_ideal)
    p_out_hot = expander_charge["state_in"].p
    h_out_hot = expander_charge["state_in"].h
    p_in_hot = p_out_hot / (1.0 - dp_recup_h)
    h_in_hot = h_out_hot + (h_out_cold - h_in_cold_actual)
    if recuperated:
        num_elements = parameters["recuperator_charge"].pop("num_elements")
    else:
        num_elements = 2
    recuperator_charge = heat_exchanger(
        working_fluid,
        h_in_hot,
        h_out_hot,
        p_in_hot,
        p_out_hot,
        working_fluid,
        h_in_cold_actual,
        h_out_cold,
        p_in_cold,
        p_out_cold,
        counter_current=True,
        num_steps=num_elements,
    )

    # Evaluate heater
    h_in_cold = expander_charge["state_out"].h
    p_in_cold = expander_charge["state_out"].p
    h_out_cold = recuperator_charge["cold_side"]["state_in"].h
    p_out_cold = recuperator_charge["cold_side"]["state_in"].p
    T_in_hot = parameters["heat_source"].pop("inlet_temperature")
    T_heat_sink_out_discharge = T_in_hot                           #an addition
    p_in_hot = parameters["heat_source"].pop("inlet_pressure")
    h_in_hot = heating_fluid.get_state(props.PT_INPUTS, p_in_hot, T_in_hot).h
    T_out_hot = heat_source_temperature_out
    p_out_hot = p_in_hot * (1 - dp_heater_h)
    h_out_hot = heating_fluid.get_state(props.PT_INPUTS, p_out_hot, T_out_hot).h
    num_elements = parameters["heater_charge"].pop("num_elements")
    heater_charge = heat_exchanger(
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
    # TODO place on other side?
    h_in = heater_charge["hot_side"]["state_out"].h
    p_in = heater_charge["hot_side"]["state_out"].p
    p_out = p_source_out
    efficiency = parameters["heat_source_pump_charge"].pop("efficiency")
    efficiency_type = parameters["heat_source_pump_charge"].pop("efficiency_type")
    heat_source_pump_charge = compression_process(
        heating_fluid,
        h_in,
        p_in,
        p_out,
        efficiency,
        efficiency_type,
    )

    # Evaluate heat sink pump
    # TODO place on other side?
    T_in = parameters["heat_sink"].pop("inlet_temperature")
    T_heat_sink_in = T_in                                         #an addition
    p_in = parameters["heat_sink"].pop("inlet_pressure")
    h_in = cooling_fluid.get_state(props.PT_INPUTS, p_in, T_in).h
    p_out = p_sink_out / (1 - dp_cooler_c)
    efficiency = parameters["heat_sink_pump_charge"].pop("efficiency")
    efficiency_type = parameters["heat_sink_pump_charge"].pop("efficiency_type")
    heat_sink_pump_charge = compression_process(
        cooling_fluid,
        h_in,
        p_in,
        p_out,
        efficiency,
        efficiency_type,
    )

    # Evaluate cooler
    p_in_cold = heat_sink_pump_charge["state_out"].p
    h_in_cold = heat_sink_pump_charge["state_out"].h
    p_out_cold = p_in_cold * (1 - dp_cooler_c)
    T_out_cold = heat_sink_temperature_out
    h_out_cold = cooling_fluid.get_state(props.PT_INPUTS, p_out_cold, T_out_cold).h
    h_in_hot = compressor_charge["state_out"].h
    p_in_hot = compressor_charge["state_out"].p
    h_out_hot = recuperator_charge["hot_side"]["state_in"].h
    p_out_hot = recuperator_charge["hot_side"]["state_in"].p
    num_elements = parameters["cooler_charge"].pop("num_elements")
    cooler_charge = heat_exchanger(
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
    W_charge = parameters.pop("charging_power") #a change 219 to 222
    m_total_charge = W_charge / (compressor_charge["specific_work"] - expander_charge["specific_work"])   #?m of working fluid?
    m_sink_charge = 1   #delete
    m_source_charge = 1   #delete
    #m*(h2-h1)=m*(h2-h1)
    #  dh_hot    dh_cold   from heat_exchanger(...)
    #m_sink_charge*(cooler["q_cold_side"]) = m_total_charge*(cooler["q_hot_side"])
    m_sink_charge = m_total_charge*(cooler_charge["q_hot_side"])/(cooler_charge["q_cold_side"])
    m_source_charge = m_total_charge*(heater_charge["q_cold_side"])/(heater_charge["q_hot_side"])
    
    # Add the mass flow to the components
    heater_charge["hot_side"]["mass_flow"] = m_source_charge
    heater_charge["cold_side"]["mass_flow"] = m_total_charge
    recuperator_charge["hot_side"]["mass_flow"] = m_total_charge
    recuperator_charge["cold_side"]["mass_flow"] = m_total_charge
    cooler_charge["hot_side"]["mass_flow"] = m_total_charge
    cooler_charge["cold_side"]["mass_flow"] = m_sink_charge
    expander_charge["mass_flow"] = m_total_charge
    compressor_charge["mass_flow"] = m_total_charge
    heat_source_pump_charge["mass_flow"] = m_source_charge
    heat_sink_pump_charge["mass_flow"] = m_sink_charge

    # Summary of components
    components = {
        "expander": expander_charge,
        "compressor": compressor_charge,
        "recuperator": recuperator_charge,
        "heater": heater_charge,
        "cooler": cooler_charge,
        "heat_source_pump": heat_source_pump_charge,
        "heat_sink_pump": heat_sink_pump_charge,
    }

    # Compute energy balances
    for name, component in components.items():

        if component["type"] == "heat_exchanger":
            hot = component["hot_side"]
            cold = component["cold_side"]
            Q1 = hot["mass_flow"] * (hot["state_in"].h - hot["state_out"].h)
            Q2 = cold["mass_flow"] * (cold["state_out"].h - cold["state_in"].h)
            component["heat_flow_rate"] = Q1
            component["heat_flow_rate_"] = Q2
            component["heat_balance"] = Q1 - Q2
            component["power"] = 0.0                                     #(want to delete)

        if component["type"] in ["compressor", "expander"]:
            component["power"] = component["mass_flow"] * component["specific_work"]
            component["heat_flow_rate"] = 0.00

    # First-law analysis
    Q_in = heater_charge["heat_flow_rate"]
    Q_out = cooler_charge["heat_flow_rate"]
    W_out = expander_charge["power"]
    W_comp = compressor_charge["power"]
    W_aux = heat_source_pump_charge["power"] + heat_sink_pump_charge["power"]
    W_in = W_comp + W_aux
    # W_net = W_out - W_in
    # Q_in_max = m_source * (heater["hot_side"]["state_in"].h - source_out_min.h)
    COP_heat_pump = Q_out / (W_in - W_out)
    COP_refrigeration = Q_in / (W_in - W_out)
    # system_efficiency = (W_out - W_in) / Q_in_max
    backwork_ratio = W_out / W_comp
    energy_balance = (Q_in + W_comp) - (W_out + Q_out)  # Ignore pumps
    W_net = compressor_charge["power"] - expander_charge["power"]

    # Define dictionary with 1st Law analysis
    energy_analysis = {
        "heater_heat_flow_charge": Q_in,
        # "heater_heat_flow_max_charge": Q_in_max,
        "recuperator_heat_flow_charge": recuperator_charge["heat_flow_rate"],
        "cooler_heat_flow_charge": Q_out,
        "expander_power_charge": W_out,
        "compressor_power_charge": W_comp,
        "heat_source_pump_power_charge": heat_source_pump_charge["power"],
        "heat_sink_pump_power_charge": heat_sink_pump_charge["power"],
        "net_cycle_power_charge": W_net,
        # "net_system_power_charge": W_out - W_in,
        "mass_flow_heating_fluid_charge": m_source_charge,
        "mass_flow_working_fluid_charge": m_total_charge,
        "mass_flow_cooling_fluid_charge": m_sink_charge,
        "COP_heat_pump": COP_heat_pump,
        "COP_refrigeration": COP_refrigeration,
        "backwork_ratio_charge": backwork_ratio,
        "energy_balance_charge": energy_balance,
    }

    # Evaluate objective function and constraints
    # output = {"components": components}
    # output = {"components": components, "energy_analysis": energy_analysis}
    # f = utilities.evaluate_objective_function(output, objective_function)
    # c_eq, c_ineq, constraint_report = utilities.evaluate_constraints(output, constraints) we have them three at the bottom
    # f = None
    # c_eq = None
    # c_ineq = None

    # Set colors for plotting
    heater_charge["hot_side"]["color"] = COLORS_MATLAB[0]
    heater_charge["cold_side"]["color"] = COLORS_MATLAB[1]
    recuperator_charge["hot_side"]["color"] = COLORS_MATLAB[1]
    recuperator_charge["cold_side"]["color"] = COLORS_MATLAB[1]
    cooler_charge["hot_side"]["color"] = COLORS_MATLAB[1]
    cooler_charge["cold_side"]["color"] = COLORS_MATLAB[6]
    expander_charge["color"] = COLORS_MATLAB[1]
    compressor_charge["color"] = COLORS_MATLAB[1]

    # Check if any fixed parameter or design variable was not used
    # utilities.check_for_unused_keys(parameters, "parameters", raise_error=True)
    # utilities.check_for_unused_keys(variables, "variables", raise_error=True)

    # Summary of components
    components = {
        "expander_charge": expander_charge,
        "compressor_charge": compressor_charge,
        "recuperator_charge": recuperator_charge,
        "heater_charge": heater_charge,
        "cooler_charge": cooler_charge,
        "heat_source_pump_charge": heat_source_pump_charge,
        "heat_sink_pump_charge": heat_sink_pump_charge,       #(want to delete)(add _charge name above before deleting)
    }

    # Cycle performance summary
    # output = {
    #     **output,
    #     "working_fluid": working_fluid,
    #     "heating_fluid": heating_fluid,
    #     "cooling_fluid": cooling_fluid,
    #     "components": components,
    #     "objective_function": f,
    #     "equality_constraints": c_eq,
    #     "inequality_constraints": c_ineq,
    # "constraints_report": constraint_report,
    # }                                          #we have at the end

    # return output   #(want to delete)

    #Written be me:
    #T_heat_source_charge = 1
    #T_heat_sink_charge = 1
    
    # T_source_in_charge = 1
    # T_source_out_charge = 1
    # T_sink_in_charge = 1
    # T_sink_out_charge = 1
    
    
    #T_heat_source_discharge = 1
    #T_heat_sink_discharge = 1
    
    T_source_in_discharge = heat_sink_temperature_out
    test = cooler_charge["cold_side"]["state_out"]["T"]
    #print(T_source_in_discharge, test)
    T_source_out_discharge = T_heat_sink_in
    T_sink_in_discharge = heat_source_temperature_out
    T_sink_out_discharge = T_heat_sink_out_discharge
    
    p_source_in_discharge = cooler_charge["cold_side"]["state_out"]["p"]
    p_source_out_discharge = cooler_charge["cold_side"]["state_in"]["p"]
    p_sink_in_discharge = heater_charge["hot_side"]["state_out"]["p"]     #think this through again
    p_sink_out_discharge = heater_charge["hot_side"]["state_in"]["p"]     #make sure it makes sense

###########################################################################################################################
#Copy and paste from cycle_power_recuperated:
# def evaluate_cycle2(
#     variables,
#     parameters,
#     constraints,
#     objective_function,
#     recuperated=True,
# ):
#     # Create copies to not change the originals
#     variables = copy.deepcopy(variables)
#     parameters = copy.deepcopy(parameters)

#     # Initialize fluid objects
#     working_fluid = props.Fluid(
#         **parameters.pop("working_fluid"), identifier="working_fluid"
#     )
#     heating_fluid = props.Fluid(
#         **parameters.pop("heating_fluid"), identifier="heating_fluid"
#     )
#     cooling_fluid = props.Fluid(
#         **parameters.pop("cooling_fluid"), identifier="cooling_fluid"   #(want to delete)
#     )

    # Extract heat source/sink parameters and give short names
    # T_source_out_max = parameters["heat_source"].pop("exit_temperature_max")
    # T_source_out_min = parameters["heat_source"].pop("exit_temperature_min")
    # T_sink_out_max = parameters["heat_sink"].pop("exit_temperature_max")
    # T_sink_out_min = parameters["heat_sink"].pop("exit_temperature_min")
    # p_source_out = parameters["heat_source"].pop("exit_pressure") #already defined
    # p_sink_out = parameters["heat_sink"].pop("exit_pressure") #already defined

    # Compute coldest state at the heat source exit
    # source_out_min = heating_fluid.get_state(
    #     props.PT_INPUTS, p_source_out, T_source_out_min
    # )

    # Extract pressure drops and give short names
    dp_heater_h = parameters["heater_discharge"].pop("pressure_drop_hot_side")
    dp_heater_c = parameters["heater_discharge"].pop("pressure_drop_cold_side")
    dp_cooler_h = parameters["cooler_discharge"].pop("pressure_drop_hot_side")
    dp_cooler_c = parameters["cooler_discharge"].pop("pressure_drop_cold_side")
    if recuperated:
        dp_recup_h = parameters["recuperator_discharge"].pop("pressure_drop_hot_side")
        dp_recup_c = parameters["recuperator_discharge"].pop("pressure_drop_cold_side")
    else:
        dp_recup_c, dp_recup_h = 0.0, 0.0

    # Extract design variables from dictionary (make sure all are used)
    expander_inlet_p = variables.pop("expander_inlet_pressure_discharge")
    expander_inlet_h = variables.pop("expander_inlet_enthalpy_discharge")
    compressor_inlet_p = variables.pop("compressor_inlet_pressure_discharge")
    compressor_inlet_h = variables.pop("compressor_inlet_enthalpy_discharge")
    # heat_source_temperature_out = variables.pop("heat_source_exit_temperature_discharge")
    # heat_sink_temperature_out = variables.pop("heat_sink_exit_temperature_discharge") #Add these ones as design variables again
    if recuperated:
        recuperator_effectiveness = variables.pop("recuperator_effectiveness_discharge")
    else:
        recuperator_effectiveness = 0.0

    # Evaluate  compressor
    dp = (1.0 - dp_heater_c) * (1.0 - dp_recup_c)
    compressor_outlet_p = expander_inlet_p / dp
    compressor_eff = parameters["compressor_discharge"].pop("efficiency")
    compressor_eff_type = parameters["compressor_discharge"].pop("efficiency_type")
    compressor_discharge = compression_process(
        working_fluid,
        compressor_inlet_h,
        compressor_inlet_p,
        compressor_outlet_p,
        compressor_eff,
        compressor_eff_type,
    )

    # Evaluate expander
    dp = (1.0 - dp_cooler_h) * (1.0 - dp_recup_h)
    expander_outlet_p = compressor_inlet_p / dp
    expander_efficiency = parameters["expander_discharge"].pop("efficiency")
    expander_efficiency_type = parameters["expander_discharge"].pop("efficiency_type")
    expander_discharge = expansion_process(
        working_fluid,
        expander_inlet_h,
        expander_inlet_p,
        expander_outlet_p,
        expander_efficiency,
        expander_efficiency_type,
    )

    # Evaluate recuperator
    eps = recuperator_effectiveness
    T_in_hot = expander_discharge["state_out"].T
    p_in_cold = compressor_discharge["state_out"].p
    h_in_cold = compressor_discharge["state_out"].h
    p_out_cold = p_in_cold * (1.0 - dp_recup_c)
    h_out_cold_ideal = working_fluid.get_state(props.PT_INPUTS, p_out_cold, T_in_hot).h
    h_out_cold_ideal = working_fluid.get_state(props.PT_INPUTS, p_out_cold, T_in_hot).h
    h_out_cold_actual = h_in_cold + eps * (h_out_cold_ideal - h_in_cold)
    p_in_hot = expander_discharge["state_out"].p
    h_in_hot = expander_discharge["state_out"].h
    p_out_hot = p_in_hot * (1.0 - dp_recup_h)
    h_out_hot = h_in_hot - (h_out_cold_actual - h_in_cold)
    if recuperated:
        num_elements = parameters["recuperator_discharge"].pop("num_elements")
    else:
        num_elements = 2
    recuperator_discharge = heat_exchanger(
        working_fluid,
        h_in_hot,
        h_out_hot,
        p_in_hot,
        p_out_hot,
        working_fluid,
        h_in_cold,
        h_out_cold_actual,
        p_in_cold,
        p_out_cold,
        counter_current=True,
        num_steps=num_elements,
    )

    # Evaluate heater
    h_in_cold = recuperator_discharge["cold_side"]["state_out"].h
    p_in_cold = recuperator_discharge["cold_side"]["state_out"].p
    h_out_cold = expander_discharge["state_in"].h
    p_out_cold = expander_discharge["state_in"].p
    T_in_hot = T_source_in_discharge   #a change
    p_in_hot = p_source_in_discharge   #a change
    h_in_hot = heating_fluid.get_state(props.PT_INPUTS, p_in_hot, T_in_hot).h
    T_out_hot = T_source_out_discharge   #a change
    p_out_hot = p_in_hot * (1 - dp_heater_h)
    h_out_hot = heating_fluid.get_state(props.PT_INPUTS, p_out_hot, T_out_hot).h
    num_elements = parameters["heater_discharge"].pop("num_elements")
    heater_discharge = heat_exchanger(
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
    # h_in = heater["hot_side"]["state_out"].h
    # p_in = heater["hot_side"]["state_out"].p
    # p_out = p_source_out
    # efficiency = parameters["heat_source_pump"].pop("efficiency")
    # efficiency_type = parameters["heat_source_pump"].pop("efficiency_type")
    # heat_source_pump_discharge = compression_process(
    #     heating_fluid,
    #     h_in,
    #     p_in,
    #     p_out,
    #     efficiency,
    #     efficiency_type,
    # )

    # Evaluate heat sink pump
    # T_in = parameters["heat_sink"].pop("inlet_temperature")
    # p_in = parameters["heat_sink"].pop("inlet_pressure")
    # h_in = cooling_fluid.get_state(props.PT_INPUTS, p_in, T_in).h
    # p_out = p_sink_out / (1 - dp_cooler_c)
    # efficiency = parameters["heat_sink_pump"].pop("efficiency")
    # efficiency_type = parameters["heat_sink_pump"].pop("efficiency_type")
    # heat_sink_pump_discharge = compression_process(
    #     cooling_fluid,
    #     h_in,
    #     p_in,
    #     p_out,
    #     efficiency,
    #     efficiency_type,
    # )

    # Evaluate cooler
    T_in_cold = T_sink_in_discharge   #an addition
    p_in_cold = p_sink_in_discharge   #an addition
    h_in_cold = cooling_fluid.get_state(props.PT_INPUTS, p_in_cold, T_in_cold).h #an addition
    # p_in_cold = heat_sink_pump["state_out"].p
    # h_in_cold = heat_sink_pump["state_out"].h
    p_out_cold = p_in_cold * (1 - dp_cooler_c)
    T_out_cold = T_sink_out_discharge   #a change
    h_out_cold = cooling_fluid.get_state(props.PT_INPUTS, p_out_cold, T_out_cold).h
    h_in_hot = recuperator_discharge["hot_side"]["state_out"].h
    p_in_hot = recuperator_discharge["hot_side"]["state_out"].p
    h_out_hot = compressor_discharge["state_in"].h
    p_out_hot = compressor_discharge["state_in"].p
    num_elements = parameters["cooler_discharge"].pop("num_elements")
    cooler_discharge = heat_exchanger(
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
    # W_net = parameters.pop("net_power")
    # expander_work = expander["specific_work"]
    # compression_work = compressor["specific_work"]
    # m_total = W_net / (expander_work - compression_work)
    # m_source = m_total * heater["mass_flow_ratio"]
    # m_sink = m_total / cooler["mass_flow_ratio"]
    m_source_discharge = m_sink_charge #an addition
    m_total_discharge = m_source_discharge*heater_discharge["q_hot_side"]/heater_discharge["q_cold_side"] #an addition
    m_sink_discharge = m_total_discharge*(cooler_discharge["q_hot_side"])/(cooler_discharge["q_cold_side"]) #an addition
    #print(m_sink_discharge, m_source_charge)
    m_error = (m_sink_discharge - m_source_charge)/m_sink_discharge

    # Add the mass flow to the components
    heater_discharge["hot_side"]["mass_flow"] = m_source_discharge
    heater_discharge["cold_side"]["mass_flow"] = m_total_discharge
    recuperator_discharge["hot_side"]["mass_flow"] = m_total_discharge
    recuperator_discharge["cold_side"]["mass_flow"] = m_total_discharge
    cooler_discharge["hot_side"]["mass_flow"] = m_total_discharge
    cooler_discharge["cold_side"]["mass_flow"] = m_sink_discharge
    expander_discharge["mass_flow"] = m_total_discharge
    compressor_discharge["mass_flow"] = m_total_discharge
    # heat_source_pump_discharge["mass_flow"] = m_source_discharge
    # heat_sink_pump_discharge["mass_flow"] = m_sink_discharge

    # Summary of components
    components = {
        **components,
        "expander_discharge": expander_discharge,
        "compressor_discharge": compressor_discharge,
        "recuperator_discharge": recuperator_discharge,
        "heater_discharge": heater_discharge,
        "cooler_discharge": cooler_discharge,
        # "heat_source_pump_discharge": heat_source_pump,
        # "heat_sink_pump_discharge": heat_sink_pump,
    }

    # Compute energy balances
    for name, component in components.items():
        if component["type"] == "heat_exchanger":
            hot = component["hot_side"]
            cold = component["cold_side"]
            Q1 = hot["mass_flow"] * (hot["state_in"].h - hot["state_out"].h)
            Q2 = cold["mass_flow"] * (cold["state_out"].h - cold["state_in"].h)
            component["heat_flow_rate"] = Q1
            component["heat_flow_rate_"] = Q2
            component["heat_balance"] = Q1 - Q2
            component["power"] = 0.0

        if component["type"] in ["compressor", "expander"]:
            component["power"] = component["mass_flow"] * component["specific_work"]
            component["heat_flow_rate"] = 0.00

    # First-law analysis
    Q_in = heater_discharge["heat_flow_rate"]
    Q_out = cooler_discharge["heat_flow_rate"]
    W_out = expander_discharge["power"]
    W_comp = compressor_discharge["power"]
    # W_aux = heat_source_pump_discharge["power"] + heat_sink_pump_discharge["power"]
    W_in = W_comp + W_aux
    # Q_in_max = m_source * (heater["hot_side"]["state_in"].h - source_out_min.h)
    cycle_efficiency = (W_out - W_in) / Q_in
    # system_efficiency = (W_out - W_in) / Q_in_max
    backwork_ratio = W_comp / W_out
    energy_balance = (Q_in + W_comp) - (W_out + Q_out)  # Ignore pumps       #commented out for now
    W_net = expander_discharge["power"] - compressor_discharge["power"]

    # Define dictionary with 1st Law analysis
    energy_analysis = {
        **energy_analysis,
        "heater_heat_flow_discharge": Q_in,
        # "heater_heat_flow_max_discharge": Q_in_max,
        "recuperator_heat_flow_discharge": recuperator_discharge["heat_flow_rate"],
        "cooler_heat_flow_discharge": Q_out,
        "expander_power_discharge": W_out,
        "compressor_power_discharge": W_comp,
        # "heat_source_pump_power_discharge": heat_source_pump_discharge["power"],
        # "heat_sink_pump_power_discharge": heat_sink_pump_discharge["power"],
        "net_cycle_power_discharge": W_net,
        "net_system_power_discharge": W_out - W_in,
        # "mass_flow_heating_fluid_discharge": m_source,
        # "mass_flow_working_fluid_discharge": m_total,
        "mass_flow_cooling_fluid_discharge": m_sink_discharge,
        "cycle_efficiency_discharge": cycle_efficiency,
        # "system_efficiency_discharge": system_efficiency,
        "backwork_ratio_discharge": backwork_ratio,
        "energy_balance_discharge": energy_balance,
    }                                                                        #commented out for now

    #RTE = W_out / W_comp#1
    RTE = energy_analysis["net_cycle_power_discharge"]/energy_analysis["net_cycle_power_charge"]
    energy_analysis["RTE"] = RTE
    #print(energy_analysis["RTE"])

    # Evaluate objective function and constraints
    output = {"components": components, "energy_analysis": energy_analysis}
    f = utilities.evaluate_objective_function(output, objective_function)
    c_eq, c_ineq, constraint_report = utilities.evaluate_constraints(output, constraints)
    c_eq.append(m_error)

    # Set colors for plotting
    heater_discharge["hot_side"]["color"] = COLORS_MATLAB[5]
    heater_discharge["cold_side"]["color"] = COLORS_MATLAB[3]
    recuperator_discharge["hot_side"]["color"] = COLORS_MATLAB[3]
    recuperator_discharge["cold_side"]["color"] = COLORS_MATLAB[3]
    cooler_discharge["hot_side"]["color"] = COLORS_MATLAB[3]
    cooler_discharge["cold_side"]["color"] = COLORS_MATLAB[2]
    expander_discharge["color"] = COLORS_MATLAB[3]
    compressor_discharge["color"] = COLORS_MATLAB[3]              #changed colours

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


