"""Cached jaxprop evaluator of the existing simple-cycle equations.

Only properties needed by the objective and constraints are calculated here.
The native model reevaluates every selected solver result before acceptance.
Uniform heater/condenser grids match case_template.yaml; no core code is changed.
"""

from thermopt.utilities.fluid_utils import get_property_fluid
import jaxprop as props

from functools import lru_cache

import run_sensitivity as run

np, th = run.np, run.th


@lru_cache(maxsize=None)
def property_fluid(name):
    return get_property_fluid(name)


def state(name, pair, a, b):
    wf = property_fluid(name)
    st = wf.get_state(pair, a, b)
    return dict(p=float(st.p), h=float(st.h), s=float(st.s), T=float(st.T))


@lru_cache(maxsize=4096)
def water_pt(p, temperature):
    return state("Water", props.PT_INPUTS, p, temperature)


@lru_cache(maxsize=4096)
def water_temperatures(h_in, h_out, p_in, p_out, count):
    return get_property_fluid('Water').get_state(props.HmassP_INPUTS, np.linspace(h_in, h_out, count), np.linspace(p_in, p_out, count)).T


def saturation(name, pressure, quality):
    wf = property_fluid(name)
    st = wf.get_state(props.PQ_INPUTS, pressure, quality)
    return dict(h=float(st.h), T=float(st.T), cp=float(st.cp))


class FastProblem(th.ThermodynamicCycleProblem):
    def fitness(self, x):
        variables = dict(zip(self.variable_names, map(float, x)))
        for key, value in variables.items():
            self.configuration["design_variables"][key]["value"] = value
        fixed = self.fixed_parameters
        name = fixed["working_fluid"]["name"]
        wf = property_fluid(name)
        pc = float(wf.critical_point.p)
        heater, cooler = fixed["heater"], fixed["cooler"]
        pump_p, pump_h = variables["compressor_inlet_pressure"], variables["compressor_inlet_enthalpy"]
        high_p, high_h = variables["expander_inlet_pressure"], variables["expander_inlet_enthalpy"]
        out_p = pump_p/(1-cooler["pressure_drop_hot_side"])
        pump_out_p = high_p/(1-heater["pressure_drop_cold_side"])
        pump_in = state(name, props.HmassP_INPUTS, pump_h, pump_p)
        pump_out_is = state(name, props.PSmass_INPUTS, pump_out_p, pump_in["s"])
        pump_work = (pump_out_is["h"]-pump_h)/fixed["compressor"]["efficiency"]
        pump_out = state(name, props.HmassP_INPUTS, pump_h+pump_work, pump_out_p)
        inlet = state(name, props.HmassP_INPUTS, high_h, high_p)
        outlet_is = state(name, props.PSmass_INPUTS, out_p, inlet["s"])
        turbine_eta = fixed["expander"]["efficiency"]
        turbine_work = turbine_eta*(high_h-outlet_is["h"])
        outlet = state(name, props.HmassP_INPUTS, high_h-turbine_work, out_p)
        sat = saturation(name, pump_p, 0)
        pump_in["subcooling"] = ((sat["h"]-pump_h)/sat["cp"] if pump_h > sat["h"] else sat["T"]-pump_in["T"])
        if name == "Propane":
            sat = saturation(name, out_p, 1)
            outlet["superheating"] = ((outlet["h"]-sat["h"])/sat["cp"] if outlet["h"] < sat["h"] else outlet["T"]-sat["T"])
        source = fixed["heat_source"]
        sink = fixed["heat_sink"]
        hot_in = water_pt(source["inlet_pressure"], source["inlet_temperature"])
        hot_out_p = source["inlet_pressure"]*(1-heater["pressure_drop_hot_side"])
        hot_out = water_pt(hot_out_p, variables["heat_source_exit_temperature"])
        source_min = water_pt(source["exit_pressure"], source["minimum_temperature"])
        source_pump_is = state("Water", props.PSmass_INPUTS, source["exit_pressure"], hot_out["s"])
        source_pump_work = (source_pump_is["h"]-hot_out["h"])/fixed["heat_source_pump"]["efficiency"]
        cold_in_p = sink["exit_pressure"]/(1-cooler["pressure_drop_cold_side"])
        cold_in_unpumped = water_pt(sink["inlet_pressure"], sink["inlet_temperature"])
        sink_pump_is = state("Water", props.PSmass_INPUTS, cold_in_p, cold_in_unpumped["s"])
        sink_pump_work = (sink_pump_is["h"]-cold_in_unpumped["h"])/fixed["heat_sink_pump"]["efficiency"]
        cold_in_h = cold_in_unpumped["h"]+sink_pump_work
        cold_out = water_pt(sink["exit_pressure"], variables["heat_sink_exit_temperature"])
        q_in, q_out = high_h-pump_out["h"], outlet["h"]-pump_h
        mass = fixed["net_power"]/(turbine_work-pump_work)
        source_mass = mass*q_in/(hot_in["h"]-hot_out["h"])
        sink_mass = mass*q_out/(cold_out["h"]-cold_in_h)
        net = fixed["net_power"]-source_mass*source_pump_work-sink_mass*sink_pump_work
        heater_n, cooler_n = heater["num_elements"], cooler["num_elements"]
        heat_hot_t = water_temperatures(hot_out["h"], hot_in["h"], hot_out_p, hot_in["p"], heater_n)
        heat_cold_t = get_property_fluid(name).get_state(props.HmassP_INPUTS, np.linspace(pump_out['h'], high_h, heater_n), np.linspace(pump_out_p, high_p, heater_n)).T
        cool_hot_t = get_property_fluid(name).get_state(props.HmassP_INPUTS, np.linspace(pump_h, outlet['h'], cooler_n), np.linspace(pump_p, out_p, cooler_n)).T
        cool_cold_t = water_temperatures(cold_in_h, cold_out["h"], cold_in_p, cold_out["p"], cooler_n)
        components = dict(
            compressor=dict(state_in=pump_in, state_out=pump_out, specific_work=pump_work),
            expander=dict(state_in=inlet, state_out=outlet, specific_work=turbine_work, pressure_ratio=high_p/out_p),
            heater=dict(temperature_difference=heat_hot_t-heat_cold_t,
                        cold_side=dict(state_in=pump_out, state_out=inlet)),
            cooler=dict(temperature_difference=cool_hot_t-cool_cold_t,
                        hot_side=dict(state_in=outlet, state_out=pump_in)))
        self.cycle_data = dict(components=components, energy_analysis=dict(
            system_efficiency=net/source_mass/(hot_in["h"]-source_min["h"])))
        self.f = th.evaluate_objective_function(self.cycle_data, self.objective_function)
        self.c_eq, self.c_ineq, self.constraint_data = th.evaluate_constraints(self.cycle_data, self.constraints)
        if name == "Propane":
            pressure = np.geomspace(min(high_p, 0.999999*pc), out_p, run.PATH_POINTS)
            ideal_h = get_property_fluid(name).get_state(props.PSmass_INPUTS, pressure, inlet['s']).h
            path_h = high_h-turbine_eta*(high_h-ideal_h)
            vapor_h = get_property_fluid(name).get_state(props.PQ_INPUTS, pressure, 1).h
            self.c_ineq = np.r_[self.c_ineq, -(path_h-vapor_h)/1e5]
        return run.psv.combine_objective_and_constraints(self.f, self.c_eq, self.c_ineq)
