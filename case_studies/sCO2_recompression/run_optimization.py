import os
import thermopt as th
import matplotlib.pyplot as plt
import jax

# Use full precision for native bicubic states and finite-difference derivatives.
jax.config.update("jax_enable_x64", True)

# Print package info
th.print_package_info()

# Define configuration filename
CONFIG_FILE = os.getenv("THERMOPT_CONFIG", "./case_sCO2_recompression.yaml")

# Initialize cycle problem
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir=os.getenv("THERMOPT_RESULTS_DIR"))
if os.getenv("DISABLE_PLOTS", "0") != "1":
    cycle.problem.plot_cycle_realtime(CONFIG_FILE)
    # Apply edited solver settings and constraint dimensions before solving.
    cycle.read_config(CONFIG_FILE)

# Perform cycle optimization
cycle.run_optimization()
cycle.save_results()

# Create an animation of the optimization progress
# cycle.create_animation(format="mp4", fps=1.0)

# Keep plots open
if os.getenv("DISABLE_PLOTS", "0") != "1":
    plt.show()
