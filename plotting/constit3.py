# Introduction
""" read constit3.tex
lambda = 1 """

# Plotting

import matplotlib.pyplot as plt
import numpy as np

alpha = [0.3, 0.6, 0.9]
colors = ["red", "blue", "green"]

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12
for spine in ["bottom", "left", "top", "right"]:
    plt.gca().spines[spine].set_zorder(0)
    plt.gca().spines[spine].set_linewidth(0.5)
    plt.gca().spines[spine].set_alpha(1)

plt.title("Constitutive Relation of Fibers with Plasticity (unbroken)", pad=13)
plt.xlabel(r"$\epsilon$", labelpad=-1)
plt.ylabel(r"$\sigma$", labelpad=8)
plt.xlim(0, 1)
plt.tick_params(axis='both', which='major', labelsize=14)

for a, c in zip(alpha, colors):
    x = np.linspace(0, 1, 1000)
    y = a*x-(1-a)*(1-np.exp(-x))
    plt.plot(x, y, color=c, linestyle='-', linewidth=1.2, zorder=10, alpha=1)

plt.legend([r"$\alpha=0.3$", r"$\alpha=0.6$", r"$\alpha=0.9$"], loc="best")
plt.savefig("plots/constit3.png", dpi=350)
