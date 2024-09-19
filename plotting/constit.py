import matplotlib.pyplot as plt  
import os  # for making output file's name the same as python file's name

# For using LaTex font in plot
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

# Read data from input file
data = []
with open("data/constit.txt", "r") as file:
    for line in file:
        x, y = line.strip().split()
        data.append((float(x), float(y)))

# Get values of x and y
x_values = [x for x, _ in data]
y_values = [y for _, y in data]

# Create plot
plt.plot(x_values, y_values, color="red", marker="none", linestyle='-', linewidth=1.2, zorder=10, alpha=1)
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$\sigma$")
plt.title("Constitutive Relation")
plt.xlim(left=0)
plt.ylim(bottom=0)

# Make plot more visible than spines in case plot is symtotically close to x axis
for spine in ["bottom", "left", "top", "right"]:
	plt.gca().spines[spine].set_zorder(0)	# lower z-order than plot
	plt.gca().spines[spine].set_linewidth(0.5)	# thinner than plot
	plt.gca().spines[spine].set_alpha(0.5)	# more transparent than plot

# Get the current script name and replace the extension with .png
script_name = os.path.basename(__file__)
plot_name = "plots/" + os.path.splitext(script_name)[0] + ".png"

# Save plot as a PNG file
plt.savefig(plot_name, dpi=350)