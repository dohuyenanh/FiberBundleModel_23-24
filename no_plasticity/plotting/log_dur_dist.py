import matplotlib.pyplot as plt 
import numpy as np 

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
for spine in ["bottom", "left", "top", "right"]:
	plt.gca().spines[spine].set_zorder(1)
	plt.gca().spines[spine].set_alpha(0.5)
	plt.gca().spines[spine].set_linewidth(0.5)

data = np.loadtxt("data/log_dur_dist.txt")
x = data[:, 0]
y = data[:, 1]

x_values = np.linspace(1, 200, 1000)

plt.title("Distribution of Profile's Duration (logarithmic binning)")
plt.plot(x, y, color="red", linestyle='-', linewidth=1.2)
plt.plot(x_values, (4*(10**7))*(x_values**(-3.95)), color="blue")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Profile's Duration = $w$")
plt.ylabel(r"$P(w)$")
plt.xlim(left=1)
plt.ylim(bottom=1)
plt.legend(["data", r"$y = (4 \cdot 10^7) \times x^{-3.95}$"], loc="best")
plt.savefig("plots/log_dur_dist.png", dpi=350)
