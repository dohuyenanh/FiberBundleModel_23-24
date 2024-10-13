import matplotlib.pyplot as plt
import json
import os  # for making output file's name the same as python file's name

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

# Read configuration
# with open("config.json", "r") as file:
#     config = json.load(file)

# Function to read data from a file
def read_data(file_path):
    data = []
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            x, y = line.strip().split()
            data.append((float(x), float(y)))
    return data

# List of input files and corresponding labels
input_files = [
    ("../data/constit(de_b-1).txt", r"$\Delta {\epsilon_b}_1=1.1$"),
    ("../data/constit(de_b-2).txt", r"$\Delta {\epsilon_b}_2=1.5$"),
    ("../data/constit(de_b-3).txt", r"$\Delta {\epsilon_b}_3=2.0$"),
    ("../data/constit(de_b-4).txt", r"$\Delta {\epsilon_b}_4=2.5$")
]

# Plot the data
for file_path, label in input_files:
    data = read_data(file_path)
    x_values, y_values = zip(*data)
    plt.plot(x_values, y_values, label=label)

# Add labels and title
plt.xlabel(r"$\epsilon$")
plt.ylabel(r"$\sigma$")
plt.title(r"Constitutive Relation with Various $\Delta \epsilon_b$")
plt.legend()
plt.xlim(0)
plt.ylim(0)

# Get the current script name and replace the extension with .png
script_name = os.path.basename(__file__)
plot_name = "plots/" + os.path.splitext(script_name)[0] + ".png"

# Save plot as a PNG file
plt.savefig(plot_name, dpi=350)