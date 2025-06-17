import numpy as np
import thermopt as th
import matplotlib.pyplot as plt

# Print package info
th.print_package_info()

# Initialize cycle problem
CONFIG_FILE = "./case_PTES_CO2.yaml"
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE)
# cycle.problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1, write_report=True)

# Perform cycle optimization
# cycle.solver.print_optimization_report(cycle.problem.x0)
cycle.run_optimization()
cycle.save_results()

# # Create an animation of the optimization progress
# #cycle.create_animation(format="mp4", fps=1)

# Keep plots open
plt.show()



