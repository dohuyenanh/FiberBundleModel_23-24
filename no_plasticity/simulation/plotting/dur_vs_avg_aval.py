import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

data = np.loadtxt("data/dur_vs_avg_aval.txt")
x = data[:, 0]
y = data[:, 1]

plt.plot(x, y, color="red", linestyle='-', linewidth=1.2, marker='o', markerfacecolor="none", markersize=5, zorder=10, alpha=1)
plt.title("Average of Avalanches Corresponding to Duration")
plt.xlabel(r"Duration = $w$")
plt.ylabel(r"Average of Avalanches = $<\Delta>$")
plt.xscale("log")
plt.yscale("log")
plt.xlim(left=1)

# Generate x values
x_values = np.linspace(1, 150, 1000)  # 1000 points from 1 to 150

plt.plot(x_values, 0.4*x_values**1.75, color="blue", alpha=0.8)
plt.legend(["data", r"$y=0.4 \times x^{1.75}$"], loc="best")

# Change the length of the x-ticks
plt.tick_params(axis='x', length=6)  

plt.savefig("plots/dur_vs_avg_aval.png", dpi=350)