import matplotlib.pyplot as plt
import numpy as np 

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
for spine in ["bottom", "left", "top", "right"]:
    plt.gca().spines[spine].set_zorder(0)
    plt.gca().spines[spine].set_linewidth(0.5)
    plt.gca().spines[spine].set_alpha(0.5)

##################################CUSTOMIZE##################################

data = np.loadtxt("data/log_aval_size.txt")
x = data[:, 0]
y = data[:, 1]

x_values = np.linspace(1, 2000, 10000)	# 10000 points from 1 to 2000

plt.plot(x, y, color="red", linestyle='-', linewidth=1.2, zorder=10, alpha=1)
plt.plot(x_values, 4*(10**6)*(x_values**(-2.5)), color="blue")
plt.title("Avalanche Distribution (logarithmic binning)")
plt.xlabel(r"Avalanche's Size = $\Delta$")
plt.ylabel(r"$P(\Delta)$")
plt.xscale("log")
plt.yscale("log")
plt.legend(["data", r"$y = (4 \cdot 10^6) \times x^{-2.5}$"], loc="best")
plt.savefig("plots/log_aval_size.png", dpi=350)
