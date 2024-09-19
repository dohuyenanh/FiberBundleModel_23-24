import matplotlib.pyplot as plt
import os

# Enable LaTeX rendering in Matplotlib
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# Define constants
a = 0.05
w_values = [17, 15, 13, 11, 9]
g = 1.03

# Read data from file
data = []
with open("data/average_profile.txt", "r") as file:
    for line in file:
        values = list(map(float, line.strip().split()))
        data.append(values)

# Transpose the data to separate columns
data = list(zip(*data))

# Plot each set of x-y values with transformations
# colors = ['#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311']
colors = ["tab:orange", "tab:purple", "tab:blue", "tab:green", "tab:red"]
labels = ['W5', 'W4', 'W3', 'W2', 'W1']

for i in range(5):
    x_values = data[2 * i]
    y_values = data[2 * i + 1]
    w = w_values[i]
    transformed_x = [x / w for x in x_values]
    transformed_y = [(y - 1) / (w ** g) for y in y_values]
    plt.plot(transformed_x, transformed_y, label=labels[i], color=colors[i], marker='o', markerfacecolor='none', linestyle='-', linewidth=1, markersize=7, alpha=0.7)

# Add labels and title with LaTeX font
plt.xlabel(r"${u}/{w}$")
plt.ylabel(r"${\langle \Delta s \rangle - 1}/{w^{\gamma}}$")
plt.title("Master Plot")

# Add legend with LaTeX font
plt.legend(loc='best')

# Set x-axis limit
plt.xlim(0, 1.1)

# Get the current script name and replace the extension with .png
script_name = os.path.basename(__file__)
plot_name = 'plots/' + os.path.splitext(script_name)[0] + '.png'

# Save plot as PNG file
plt.savefig(plot_name, dpi=350)
