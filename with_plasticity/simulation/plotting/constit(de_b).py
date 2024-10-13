import matplotlib.pyplot as plt
import os  # for making output file's name the same as python file's name

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

# Function to read data from a file
def read_data(file_path):
    data = []
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            x, y = line.strip().split()
            data.append((float(x), float(y)))
    return data

# Read data from input files
data1 = read_data("../data/constit(de_b-1).txt")
data2 = read_data("../data/constit(de_b-2).txt")
data3 = read_data("../data/constit(de_b-3).txt")

# Get values of x and y for each dataset
x_values1, y_values1 = zip(*data1) 
x_values2, y_values2 = zip(*data2)
x_values3, y_values3 = zip(*data3)

# Plot the data
plt.plot(x_values1, y_values1, label=r"$\Delta {\epsilon_b}_1=1.1$", color="red")
plt.plot(x_values2, y_values2, label=r"$\Delta {\epsilon_b}_2=1.5$", color="green")
plt.plot(x_values3, y_values3, label=r"$\Delta {\epsilon_b}_3=2.0$", color="blue")

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