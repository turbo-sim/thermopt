import numpy as np
import thermopt as th
import matplotlib.pyplot as plt


#Varying isentropic efficiency
CONFIG_FILE = "./case_PTES_CO2.yaml"
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE)
#cycle.problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1, write_report=True) #maybe not needed, because we set the yaml file once
RTE = []
isentropic_efficiency = np.arange(0.7, 1.01, 0.05) #from 0.5 and 0.6 created an error in calculations
x0 = None #empty space
for i in isentropic_efficiency:
    i = round(i, 2)
    print(i)
    # Run optimization and store results
    cycle.set_config_value("problem_formulation.fixed_parameters.expander_charge.efficiency", i)
    cycle.set_config_value("problem_formulation.fixed_parameters.compressor_charge.efficiency", i)
    cycle.set_config_value("problem_formulation.fixed_parameters.expander_discharge.efficiency", i)
    cycle.set_config_value("problem_formulation.fixed_parameters.compressor_discharge.efficiency", i)
    cycle.run_optimization(x0=x0) #x0 stands for first value (initial guess)
    cycle.save_results()
    plt.close(cycle.problem.figure)
    x0 = cycle.solver.x_final
    # Extract results
    RTripE = cycle.problem.cycle_data["energy_analysis"]["roundtrip_efficiency"]
    #Store
    RTE.append(RTripE)
plt.plot(isentropic_efficiency, RTE)
plt.xlabel('machine isentropic efficiency')
plt.ylabel('round-trip efficiency')
plt.show()