import numpy as np
import thermopt as th
import matplotlib.pyplot as plt

# Pressure ratio range
Pi = np.linspace(1.01, 50, 500)

# Efficiency values to test
eta_list = [0.5, 0.8, 0.9]
gamma_list = [1.3, 1.4, 1.5]
linestyles = ["-", "--", ":"]  # One for each gamma

# Define entropy generation from isentropic efficiency
def entropy_gen_from_eta_s(Pi, eta_s, gamma):
    Pi_inv = 1 / Pi
    exponent = (gamma - 1) / gamma
    T_ratio = (1 - eta_s * (1 - Pi_inv**exponent))**(1 / exponent)
    s_gen = np.log(T_ratio / Pi_inv)
    return s_gen

# Define entropy generation from polytropic efficiency
def entropy_gen_polytropic_expander(Pi, eta_p):
    return (1 - eta_p) * np.log(Pi)

# Setup plot
plt.figure(figsize=(7, 5))
plt.xlabel("Pressure ratio $\\Pi$")
plt.ylabel("Entropy generation $s_{\\mathrm{gen}}$")
plt.grid(True)

# Color for each efficiency
colors = plt.get_cmap("magma")(np.linspace(0.2, 0.8, len(eta_list)))

for i, eta in enumerate(eta_list):
    color = colors[i]
    for j, gamma in enumerate(gamma_list):
        ls = linestyles[j]
        s_gen = entropy_gen_from_eta_s(Pi, eta, gamma)
        plt.plot(Pi, s_gen, color=color, linestyle=ls, label=f"$\\eta_s$={eta}, $\\gamma$={gamma}")

    # Polytropic curve (no gamma dependence), using dash-dot style
    s_gen_p = entropy_gen_polytropic_expander(Pi, eta)
    plt.plot(Pi, s_gen_p, color=color, linestyle="-.", linewidth=2,
             label=f"$\\eta_p$={eta} (no $\\gamma$)")

# Finalize
plt.legend(fontsize=9, ncol=1)
plt.tight_layout()
plt.show()