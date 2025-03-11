import numpy as np
import thermopt as th
import matplotlib.pyplot as plt


# Create fluid
fluid_name = "CO2"
fluid = th.Fluid(name=fluid_name, exceptions=True)

# Define inlet state and exit pressure
temperature_in = 500
pressure_in = 2 * fluid.critical_point.p
pressure_out = pressure_in/2

# Compute inlet state
state_in = fluid.get_state(th.PT_INPUTS, pressure_in, temperature_in)
entropy_in = state_in.s

# Compute properties along isentropic expansion
p_list = np.linspace(pressure_in, pressure_out, 10)
T_list = []
s_list = []
for p in p_list:
    state = fluid.get_state(th.PSmass_INPUTS, p, entropy_in)
    T_list.append(state.T)
    s_list.append(state.s)

# Plot the states in T-s coordinates
fig, ax = plt.subplots()
ax.set_xlabel("Entropy (J/kg/K)")
ax.set_ylabel("Temperature (K)")
ax.plot(s_list, T_list, "ko")

# Plot the phase diagram
fluid.plot_phase_diagram(x_prop="s", y_prop="T", axes=ax,
                         plot_critical_point=True,
                         plot_saturation_line=True,
                         plot_quality_isolines=True)

fig.tight_layout(pad=1)

plt.show()



