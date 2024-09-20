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

data = np.loadtxt("data/lin_aval_size.txt")
x = data[:, 0]
y = data[:, 1]

x_values = np.linspace(1, 1000, 10000)	# 10000 points from 1 to 1000

plt.plot(x, y, color="red", linestyle='-', linewidth=1.2, zorder=10, alpha=1)
plt.plot(x_values, (10**7)*(x_values**(-2.5)), color="blue")
plt.title("Avalanche Distribution (linear binning)")
plt.xlabel(r"Avalanche's Size = $\Delta$")
plt.ylabel(r"$P(\Delta)$")
plt.xscale("log")
plt.yscale("log")
plt.xlim(left=1)
plt.ylim(bottom=1)
plt.legend(["data", r"$y=10^7 \times x^{-2.5}$"], loc="best")
plt.savefig("plots/lin_aval_size.png", dpi=350)
