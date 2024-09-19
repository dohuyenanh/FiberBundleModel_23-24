# Introduction
""" 
Suppose the probability density functions of epsilon_1 and epsilon_2 are both continuos uniform distribution, i.e. their PDF and CDF look like:
       f(epsilon_1) = 1/(2*de1)
       F(epsilon_1) = (epsilon - e1 + de1)/(2*de1)
       g(epsilon_2) = 1/(2*de2)
       G(epsilon_2) = (epsilon - e2 + de2)/(2*de2)

Average number of intact fibers:    1 - F(E*epsilon)
Average number of yielding fibers:  F(E*epsilon) * [1 - G(E*epsilon)]
Average number of broken fibers:    F(E*epsilon) * G(E*epsilon)
"""

# Known quantities:
E = 1       # elastic modulus
e1 = 0.3    # mean value of epsilon_1
e2 = 0.6    # mean value of epsilon_2
de1 = 0.05  # distance between epsilon_1 and the 2 ends
de2 = 0.05  # distance between epsilon_2 and the 2 ends 

# Plotting

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
for spine in ["bottom", "left", "top", "right"]:
    plt.gca().spines[spine].set_zorder(0)
    plt.gca().spines[spine].set_linewidth(0.5)
    plt.gca().spines[spine].set_alpha(1)
    
plt.title("Constitutive Relation of Fibers with Plasticity", pad=13)
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$\sigma$")
plt.xlim(0, 1)
    
# intact region
x = np.linspace(0, e1, 1000)
y = 1-(x-e1+de1)/(2*de1)
plt.plot(x, y, color="red", linestyle='-', linewidth=1.2, zorder=10, alpha=1)

# yielding region
x = np.linspace(e1, e2, 1000)
y = 1-(x-e2+de2)/(2*de2)
plt.plot(x, y, color="blue", linestyle='-', linewidth=1.2, zorder=10, alpha=1)

# broken region
x = np.linspace(e2, 1, 1000)
y = 0*x
plt.plot(x, y, color="green", linestyle='-', linewidth=1.2, zorder=10, alpha=1)

plt.legend(["Intact", "Yielding", "Broken"], loc="best")
plt.savefig("plots/constit2.png", dpi=350)