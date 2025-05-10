import thermopt as th
import matplotlib.pyplot as plt


# Print package info
th.print_package_info()

# Read configuration file
CONFIG_FILE = "./case_sCO2_PTES_roberto.yaml"
# CONFIG_FILE = "./case_sCO2_heat_pump_charge.yaml"
# CONFIG_FILE = "./case_sCO2_heat_pump_air.yaml"

config = th.read_configuration_file(CONFIG_FILE)

# Initialize cycle problem
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE)
cycle.problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1, write_report=True)

# Perform cycle optimization
cycle.run_optimization()
cycle.save_results()

# Create an animation of the optimization progress
#cycle.create_animation(format="mp4", fps=1)

# Keep plots open
plt.show()
