import numpy as np
import thermopt as th
import matplotlib.pyplot as plt


# Pressure ratio range
Pi = np.linspace(1.01, 50, 500)

# Polytropic efficiencies to analyze
eta_p_list = [0.9, 0.8, 0.7]

# Define entropy generation function (same for both)
def entropy_gen_polytropic_compressor(Pi, eta_p):
    return (1 / eta_p - 1) * np.log(Pi)

# Define entropy generation function (same for both)
def entropy_gen_polytropic_turbine(Pi, eta_p):
    return (1 - eta_p) * np.log(Pi)

# Set up colors and styles
colors = plt.get_cmap("magma")(np.linspace(0.2, 0.8, len(eta_p_list)))

# Initialize plot
plt.figure(figsize=(7, 5))
plt.xlabel("Pressure ratio $\\Pi$")
plt.ylabel("Entropy generation $s_{\\mathrm{gen}} / R$")
plt.grid(True)

# Plot compressor cases
for i, eta_p in enumerate(eta_p_list):
    s_gen_comp = entropy_gen_polytropic_compressor(Pi, eta_p)
    plt.plot(Pi, s_gen_comp, color=colors[i], linestyle='-', label=f"$\\eta_{{\\mathrm{{comp}}}} = {eta_p}$")
    
# Plot turbine cases
for i, eta_p in enumerate(eta_p_list):
    s_gen_turb = entropy_gen_polytropic_turbine(Pi, eta_p)
    plt.plot(Pi, s_gen_turb, color=colors[i], linestyle='--', label=f"$\\eta_{{\\mathrm{{turb}}}} = {eta_p}$")

# Finalize
plt.legend(fontsize=11, ncol=2)
plt.tight_layout()
plt.show()
