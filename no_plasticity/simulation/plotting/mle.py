import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.ticker import LogLocator, ScalarFormatter, FuncFormatter

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 10
for spine in ["bottom", "left", "top", "right"]:
    plt.gca().spines[spine].set_zorder(0)
    plt.gca().spines[spine].set_linewidth(0.5)
    plt.gca().spines[spine].set_alpha(1)

data = np.loadtxt("data/alpha(x_min).txt")
x = data[:, 0]
y = data[:, 1]

x_values = np.linspace(1, 125, 10000)	# 10000 points from 1 to 1000

plt.plot(x, y, color="red", linestyle='-', linewidth=1.2, marker='o', markerfacecolor="none", markersize=5, zorder=10, alpha=1)
plt.plot(x_values, x_values*0+0.999999, color="blue")
plt.title("Mean of the maximum likelihood estimate for the scaling parameter", fontdict={"fontsize": 14}, pad=15)
plt.xlabel(r"$x_{min}$")
plt.ylabel(r"$\alpha$")
plt.xscale("log")
# plt.yscale("log")
plt.xlim(1, 200)
plt.ylim(0.999990, max(y)+0.000005)
plt.tick_params(axis='x', which='major', length=7)
plt.tick_params(axis='x', which='minor', length=4)
plt.tick_params(axis='y', labelsize=9)

""" plt.grid(True, which="both", axis='x', ls="--", alpha=0.5)
# Customize y-axis tick steps
y_ticks = np.arange(0.999990, max(y), 0.000005)  # ticks from [0] to [1] with step [2]
plt.yticks(y_ticks, [f"{tick:.6f}" for tick in y_ticks])  # Format ticks as plain numbers """

# Get current Axes instance
ax = plt.gca()

# Use ScalarFormatter to avoid scientific notation
ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

# Set major ticks on x-axis
ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=3))

# Set minor ticks on x-axis
ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs="all", numticks=12))
ax.xaxis.set_minor_formatter(FuncFormatter(lambda x, _: '2' if (x / 10**int(np.log10(x))) == 2 else ('5' if (x / 10**int(np.log10(x))) == 5 else '')))

plt.savefig("plots/mle.png", dpi=350)
