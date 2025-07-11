
import numpy as np
import thermopt as th
import matplotlib.pyplot as plt

# Define the isentropic efficiency for a compressor
def compressor_eta_s(Pi, eta_p, gamma):
    exponent = (gamma - 1) / gamma
    numerator = Pi**exponent - 1
    denominator = Pi**(exponent / eta_p) - 1
    return numerator / denominator

# Define the isentropic efficiency for a turbine
def turbine_eta_s(Pi, eta_p, gamma):
    Pi = 1 / Pi
    exponent = (gamma - 1) / gamma
    numerator = 1 - Pi**(eta_p * exponent)
    denominator = 1 - Pi**exponent
    return numerator / denominator


# Define parameter range
Pi = np.linspace(1.00001, 30, 500)  # Avoid Pi = 1 to prevent division by zero
eta_p_list = [0.6, 0.7, 0.8, 0.9, 1.0]
gamma_list = [1.3, 1.4, 1.5]
linestyles = ["-", "--", ":"]

# Plot compressor trends
plt.figure(figsize=(6, 5))
plt.ylim([0.2, 1.2])
plt.xlabel("Pressure ratio $\\Pi$")
plt.ylabel("Compressor isentropic efficiency $\\eta_s$")
colors = plt.get_cmap("magma")(np.linspace(0.2, 0.8, len(eta_p_list)))
for i, eta_p in enumerate(eta_p_list):
    for j, gamma in enumerate(gamma_list):
        eta_s = compressor_eta_s(Pi, eta_p, gamma)
        label = f"$\\gamma$ = {gamma}" if i == 0 else None  # Show gamma label only once
        plt.plot(Pi, eta_s, color=colors[i], linestyle=linestyles[j], label=label)
plt.legend(loc="upper right", ncol=len(gamma_list))
plt.tight_layout(pad=1)

# Plot turbine trends
plt.figure(figsize=(6, 5))
plt.ylim([0.2, 1.2])
plt.xlabel("Pressure ratio $\\Pi$")
plt.ylabel("Turbine isentropic efficiency $\\eta_s$")
colors = plt.get_cmap("magma")(np.linspace(0.2, 0.8, len(eta_p_list)))
for i, eta_p in enumerate(eta_p_list):
    for j, gamma in enumerate(gamma_list):
        eta_s = turbine_eta_s(Pi, eta_p, gamma)
        label = f"$\\gamma$ = {gamma}" if i == 0 else None  # Only once per color
        plt.plot(Pi, eta_s, color=colors[i], linestyle=linestyles[j], label=label)
plt.legend(loc="upper right", ncol=len(gamma_list))
plt.tight_layout(pad=1)


# Show figures
plt.show()

