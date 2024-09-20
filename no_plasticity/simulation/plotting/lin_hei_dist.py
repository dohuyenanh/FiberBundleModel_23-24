import matplotlib.pyplot as plt 
import numpy as np 

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
for spine in ["bottom", "left", "top", "right"]:
	plt.gca().spines[spine].set_zorder(1)
	plt.gca().spines[spine].set_alpha(0.5)
	plt.gca().spines[spine].set_linewidth(0.5)

data = np.loadtxt("data/lin_hei_dist.txt")
x = data[:, 0]
y = data[:, 1]

x_values = np.linspace(1, 70, 1000)

plt.title("Profile's Height Distribution (linear binning)")
plt.plot(x, y, color="red", linestyle='-', linewidth=1.2, marker='o', markerfacecolor="none", markersize=5, zorder=10, alpha=1)
plt.plot(x_values, 1*x_values**-4, color="blue")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Profile's Height = $h$")
plt.ylabel(r"$P(h)$")
plt.legend(["data", r"$y = x^{-4}$"], loc="best")
plt.savefig("plots/lin_hei_dist.png", dpi=350)
