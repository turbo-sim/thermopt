"""Run one simple CO2 case; the pressure-ratio sweep is a separate program."""

import os
import jax
import matplotlib.pyplot as plt
import thermopt as th

jax.config.update("jax_enable_x64", True)
CONFIG_FILE = os.getenv("THERMOPT_CONFIG", "case_sCO2_simple.yaml")
cycle = th.ThermodynamicCycleOptimization(
    CONFIG_FILE, out_dir=os.getenv("THERMOPT_RESULTS_DIR")
)
cycle.run_optimization()
cycle.save_results()
if os.getenv("DISABLE_PLOTS", "0") != "1":
    plt.show()
