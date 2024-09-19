import matplotlib.pyplot as plt
import os

# Read data from file
data = []
with open("data/wait_time_dist.txt", "r") as file:
    for line in file:
        x, y = line.strip().split()
        data.append((float(x), float(y)))

# Extract x and y values
x_values = [x for x, _ in data]
y_values = [y for _, y in data]

# Create plot with higher z-order for the data
plt.plot(x_values, y_values, color='red', marker='none',
         linestyle='-', linewidth=2, markersize=2, zorder=5, alpha=0.8)
plt.xlabel(r"$\Delta \sigma (\div 10^6)$")
plt.ylabel(r"$P(\Delta \sigma)$")
plt.title("Distribution of Waiting Time Between Avalanches")

# Adjust the z-order of the x-axis
plt.gca().spines['bottom'].set_zorder(0)

# Set the limits to ensure both axes start at 0
plt.xlim(left=0)
plt.ylim(bottom=0)

""" # Set the x-axis and y-axis to logarithmic scale
plt.xscale('log')
plt.yscale('log') """

# Set x-axis ticks with a step of 2
plt.xticks(range(0, int(max(x_values)) + 1, 2))

# Get the current script name and replace the extension with .png
script_name = os.path.basename(__file__)
plot_name = 'plots/' + os.path.splitext(script_name)[0] + '.png'

# Save plot as PNG file
plt.savefig(plot_name, dpi=350)

""" # Show the plot
plt.show() """
