import numpy as np
import thermopt as th
import matplotlib.pyplot as plt




#Varying hot storage maximum temperature for different values of discharge maximum pressure
CONFIG_FILE = "./case_PTES_CO2.yaml"
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE)
#cycle.problem.plot_cycle_realtime(CONFIG_FILE, update_interval=0.1, write_report=True) #maybe not needed, because we set the yaml file once
RTEtemp = []
hot_storage_maximum_temperature = np.arange(400 + 273.15, 801 + 273.15, 50) #from 200 created an error in calculations
discharge_maximum_pressure = [150e5, 200e5, 250e5]
RTEtempmatrix = np.zeros((len(discharge_maximum_pressure),len(hot_storage_maximum_temperature)))
for p in discharge_maximum_pressure:
    print(p)
    #where is maximum pressure?
    cycle.set_config_value("problem_formulation.design_variables.expander_inlet_pressure_charge.max", p)
    cycle.set_config_value("problem_formulation.design_variables.expander_inlet_pressure_discharge.max", p)
    x0 = None #empty space
    for i in hot_storage_maximum_temperature:
        i = round(i, 2)
        print(i)
        # Run optimization and store results
        cycle.set_config_value("problem_formulation.design_variables.hot_storage_upper_temperature.min", i)
        cycle.set_config_value("problem_formulation.design_variables.hot_storage_upper_temperature.max", i)
        cycle.set_config_value("problem_formulation.design_variables.hot_storage_upper_temperature.value", i)
        cycle.run_optimization(x0=x0) #x0 stands for first value (initial guess)
        cycle.save_results()
        plt.close(cycle.problem.figure)
        x0 = cycle.solver.x_final
        # Extract results
        RTripEtemp = cycle.problem.cycle_data["energy_analysis"]["roundtrip_efficiency"]
        #Store
        RTEtemp.append(RTripEtemp)
    RTEtempmatrix[discharge_maximum_pressure.index(p)] = RTEtemp
plt.plot(hot_storage_maximum_temperature, RTEtempmatrix[0])
plt.plot(hot_storage_maximum_temperature, RTEtempmatrix[1])
plt.plot(hot_storage_maximum_temperature, RTEtempmatrix[2])
plt.xlabel('hot storage maximum temperature')
plt.ylabel('round-trip efficiency')
plt.legend(["150bar","200bar","250bar"], loc="lower left")
plt.show()
plt.plot(hot_storage_maximum_temperature, RTEtemp)
plt.xlabel('hot storage maximum temperature')
plt.ylabel('round-trip efficiency')
plt.show()


