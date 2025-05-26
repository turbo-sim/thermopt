import thermopt as th
import matplotlib.pyplot as plt

# Print package info
# th.print_package_info()

# Read configuration file
CONFIG_FILE = "./case_PTES_CO2.yaml"
config = th.read_configuration_file(CONFIG_FILE)
problem = th.ThermodynamicCycleProblem(config["problem_formulation"])
problem.fitness(problem.x0) #fitness: the same as calling evaluate_cycle function
#print(problem.cycle_data)
problem.plot_cycle()
#problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1)
# print(problem.cycle_data["energy_analysis"])

# # Initialize cycle problem
# cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE)
# cycle.problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1, write_report=True)

# # Perform cycle optimization
# cycle.run_optimization()
# cycle.save_results()

# # Create an animation of the optimization progress
# #cycle.create_animation(format="mp4", fps=1)

# # Keep plots open
# plt.show()
