import os
import thermopt as th
import matplotlib.pyplot as plt

# Print package info
th.print_package_info()

# Initialize cycle problem
CONFIG_FILE = "./case_water.yaml"
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir=os.getenv("THERMOPT_RESULTS_DIR"))

# Only show interactive plots when not in test mode
if os.getenv("DISABLE_PLOTS", "0") != "1":
    cycle.problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1)
    # Apply edited solver settings and constraint dimensions before solving.
    cycle.read_config(CONFIG_FILE)

# Perform cycle optimization
cycle.run_optimization()
cycle.save_results()

# Create an animation of the optimization progress
if os.getenv("DISABLE_PLOTS", "0") != "1":
    cycle.create_animation(format="mp4", fps=1)

# Show figures
if os.getenv("DISABLE_PLOTS", "0") != "1":
    plt.show()


