import os
import thermopt as th
import matplotlib.pyplot as plt

# Print package info
th.print_package_info()

# Define configuration filename
config_files = ["./case_butane_ORC.yaml"]
# config_files = ["./case_butane_PEORC.yaml"]
# config_files = ["./case_butane_transcritical.yaml"]
# config_files = ["./case_butane_ORC.yaml", "./case_butane_PEORC.yaml","./case_butane_transcritical.yaml"]
# config_files = ["./case_Novec649_PEORC.yaml"]

# Loop over all cases
for config in config_files:

    # Initialize Brayton cycle problem
    cycle = th.ThermodynamicCycleOptimization(config, out_dir=os.getenv("THERMOPT_RESULTS_DIR"))
    if os.getenv("DISABLE_PLOTS", "0") != "1":
        cycle.problem.plot_cycle_realtime(config, update_interval=0.1)

    # Perform cycle optimization
    cycle.run_optimization()
    cycle.save_results()

    # # Create an animation of the optimization progress
    # cycle.create_animation(format="mp4", fps=1)

# Keep plots open
if os.getenv("DISABLE_PLOTS", "0") != "1":
    plt.show()


