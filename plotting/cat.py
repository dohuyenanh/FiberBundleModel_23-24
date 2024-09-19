import matplotlib.pyplot as plt
import os
import numpy as np

# Enable LaTeX rendering in Matplotlib
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Read data from file
data = np.loadtxt("data/dist_cat_aval.txt")

# Extract x and y values
x_values = data[:, 0]   # ':' means all rows, '0' means first column
y_values = data[:, 1]   # ':' means all rows, '1' means second column

# Create plot with higher z-order for the data
plt.plot(x_values, y_values, color='red', marker='none',
         linestyle='-', linewidth=1.2, markersize=2, zorder=5, alpha=1)
plt.xlabel("Catastrophic Avalanche's Size = x")
plt.ylabel("P(x)")
plt.title("Catastrophic Avalanche Distribution")

# Adjust the z-order, linewidth, and alpha for all spines in case plot is symtotically close to the axes
for spine in ['bottom', 'left', 'top', 'right']:
    plt.gca().spines[spine].set_zorder(0)
    plt.gca().spines[spine].set_linewidth(0.5)
    plt.gca().spines[spine].set_alpha(0.5)

# Set the limits to ensure both axes start at 0 (or not)
plt.xlim(right=50000)
plt.ylim(bottom=0, top=0.0008)

# Get the current script name and replace the extension with .png
script_name = os.path.basename(__file__)
plot_name = 'plots/' + os.path.splitext(script_name)[0] + '.png'

# Save plot as PNG file
plt.savefig(plot_name, dpi=350)
