import os
import jax
import thermopt as th
import matplotlib.pyplot as plt

# Native bicubic states and finite-difference derivatives require full precision.
jax.config.update("jax_enable_x64", True)

# Print package info
th.print_package_info()

# Define configuration filename
CONFIG_FILE = os.getenv("THERMOPT_CONFIG", "case_sCO2_recuperated.yaml")

# Initialize cycle problem
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir=os.getenv("THERMOPT_RESULTS_DIR"))
if os.getenv("DISABLE_PLOTS", "0") != "1":
    cycle.problem.plot_cycle_realtime(CONFIG_FILE)
    # Apply edited solver settings and constraint dimensions before solving.
    cycle.read_config(CONFIG_FILE)

# Optimize cycle
cycle.run_optimization()
# cycle.save_results()

# x_final_norm = cycle.solver.x_final_norm
# cycle.solver.make_constraint_report(x_final_norm, 1e-5)
cycle.solver.print_optimization_report(
    tol=1e-6,
    include_kkt_conditions=True,
    include_design_variables=True,
    include_constraints=True,
    include_multipliers=True,
)

# # Create an animation of the optimization progress
# cycle.create_animation(format="mp4", fps=1.0)

# Keep plots open
if os.getenv("DISABLE_PLOTS", "0") != "1":
    plt.show()
