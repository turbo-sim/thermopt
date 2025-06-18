import numpy as np
import matplotlib.pyplot as plt
import thermopt as th
th.print_package_info()

# -------------------------------------------------------------------------- #
# First run: Optimization with turbomachinery models (geometric constraints)
# -------------------------------------------------------------------------- #
CONFIG_FILE = "../case_PTES_CO2_turbo.yaml"
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir="results_turbo_model")

# Optional: live plotting of the cycle during optimization
# cycle.problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1, write_report=True)
cycle.set_config_value("solver_options.callbacks.plot_cycle", True)

# Run optimization
cycle.run_optimization()
cycle.save_results()

# -------------------------------------------------------------------------- #
# Extract final isentropic efficiencies from turbo model
# -------------------------------------------------------------------------- #
eta_compressor_charge = cycle.problem.cycle_data["components"]["compressor_charge"]["data_out"]["isentropic_efficiency"]
eta_compressor_discharge = cycle.problem.cycle_data["components"]["compressor_discharge"]["data_out"]["isentropic_efficiency"]
eta_expander_charge = cycle.problem.cycle_data["components"]["expander_charge"]["data_out"]["isentropic_efficiency"]
eta_expander_discharge = cycle.problem.cycle_data["components"]["expander_discharge"]["data_out"]["isentropic_efficiency"]

print("\nFinal isentropic efficiencies from turbo model:")
print(f"  Compressor charge     eta_is = {eta_compressor_charge:.6f}")
print(f"  Compressor discharge  eta_is = {eta_compressor_discharge:.6f}")
print(f"  Expander charge       eta_is = {eta_expander_charge:.6f}")
print(f"  Expander discharge    eta_is = {eta_expander_discharge:.6f}")

# -------------------------------------------------------------------------- #
# Second run: Fixed-efficiency model (no turbomachinery constraints)
# -------------------------------------------------------------------------- #
CONFIG_FILE = "../case_PTES_CO2.yaml"
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir="results_fixed_efficiency")
cycle.set_config_value("solver_options.callbacks.plot_cycle", True)

# Apply the final efficiencies from turbo model
cycle.set_config_value("problem_formulation.fixed_parameters.compressor_charge.efficiency", eta_compressor_charge)
cycle.set_config_value("problem_formulation.fixed_parameters.compressor_discharge.efficiency", eta_compressor_discharge)
cycle.set_config_value("problem_formulation.fixed_parameters.expander_charge.efficiency", eta_expander_charge)
cycle.set_config_value("problem_formulation.fixed_parameters.expander_discharge.efficiency", eta_expander_discharge)

# Run and save
cycle.run_optimization()
cycle.save_results()





# # -------------------------------------------------------------------------- #
# # -------------------------------------------------------------------------- #
# CONFIG_FILE = "../case_PTES_CO2_turbo.yaml"
# cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir="output_constrained")
# cycle.set_config_value("solver_options.callbacks.plot_cycle", True)


# # Set temperature and pressure constraints
# T_hot = 600 + 273.15
# var_temp = "problem_formulation.design_variables.hot_storage_upper_temperature"
# cycle.set_config_value(f"{var_temp}.max", T_hot)
# cycle.set_config_value(f"{var_temp}.value", T_hot)




# cycle.set_constraint(
#     variable="$components.expander_charge.data_out.angular_speed - $components.compressor_charge.data_out.angular_speed",
#     type="=",
#     value=0.0,
#     normalize=1000,
# )
# cycle.set_constraint(
#     variable="$components.expander_discharge.data_out.angular_speed - $components.compressor_discharge.data_out.angular_speed",
#     type="=",
#     value=0.0,
#     normalize=1000,
# )

# cycle.run_optimization(x0=x_final)
# cycle.save_results()




# print()
# print("Expander charge")
# th.print_dict(cycle.problem.cycle_data["components"]["expander_charge"]["data_in"])
# th.print_dict(cycle.problem.cycle_data["components"]["expander_charge"]["data_out"])

# print()
# print("Compressor charge")
# th.print_dict(cycle.problem.cycle_data["components"]["compressor_charge"]["data_in"])
# th.print_dict(cycle.problem.cycle_data["components"]["compressor_charge"]["data_out"])

# print()
# print("Expander discharge")
# th.print_dict(cycle.problem.cycle_data["components"]["expander_discharge"]["data_in"])
# th.print_dict(cycle.problem.cycle_data["components"]["expander_discharge"]["data_out"])

# print()
# print("Compressor discharge")
# th.print_dict(cycle.problem.cycle_data["components"]["compressor_discharge"]["data_in"])
# th.print_dict(cycle.problem.cycle_data["components"]["compressor_discharge"]["data_out"])



# # cycle.save_results()

# # # Create an animation of the optimization progress
# # #cycle.create_animation(format="mp4", fps=1)

# # Keep plots open
# # plt.show()W


# # # Set temperature and pressure constraints (unreliable convergence?)
# # T_hot = 600 + 273.15
# # var_temp = "problem_formulation.design_variables.hot_storage_upper_temperature"
# # cycle.set_config_value(f"{var_temp}.max", T_hot)
# # cycle.set_config_value(f"{var_temp}.value", T_hot)





# # # Hardcoded turbomachinery parameters
# # specific_speed = 0.7         # nondimensional
# # flow_coefficient = 0.06      # nondimensional
# # backsweep_angle = 25.0       # degrees

# # # Define variable names and set values
# # variables_to_set = {
# #     # Charge side
# #     "compressor_flow_coefficient_charge": flow_coefficient,
# #     "compressor_backsweep_charge": backsweep_angle,
# #     "expander_specific_speed_charge": specific_speed,

# #     # Discharge side
# #     "compressor_flow_coefficient_discharge": flow_coefficient,
# #     "compressor_backsweep_discharge": backsweep_angle,
# #     "expander_specific_speed_discharge": specific_speed,
# # }

# # # Set min, max, and value to the same number
# # for var_name, value in variables_to_set.items():
# #     base_key = f"problem_formulation.design_variables.{var_name}"
# #     cycle.set_config_value(f"{base_key}.min", value)
# #     cycle.set_config_value(f"{base_key}.max", value)
# #     cycle.set_config_value(f"{base_key}.value", value)

