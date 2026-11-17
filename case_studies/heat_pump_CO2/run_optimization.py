import os
import thermopt as th
import matplotlib.pyplot as plt

# Print package info
th.print_package_info()

# Define configuration filename
CONFIG_FILE = "./case_sCO2_heat_pump_valve.yaml"
# CONFIG_FILE = "./case_sCO2_heat_pump_turbine.yaml"

# Initialize cycle problem
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir=os.getenv("THERMOPT_RESULTS_DIR"))
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

# Keep plots open
if os.getenv("DISABLE_PLOTS", "0") != "1":
    plt.show()


