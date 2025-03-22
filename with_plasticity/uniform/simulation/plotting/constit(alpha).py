import os
import matplotlib.pyplot as plt
from read_parameters import read_parameters, read_number_of_parameters  # pylint: disable=import-error

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

# Read alpha parameters and number of alpha values from the JSON file
alpha_params = read_parameters("../program/parameters.json", "alpha_1", "alpha_2", "alpha_3")
no_of_alpha = read_number_of_parameters("../program/parameters.json", "no. of alpha")

# List of input files and corresponding labels
input_files = [
    (f"../data/constit(alpha-{i+1}).txt", rf"$\alpha_{i+1}={alpha_params[f'alpha_{i+1}']}$")
    for i in range(no_of_alpha)
]

# Function to read data from a file
def read_data(file_path):
    data = []
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            x, y = line.strip().split()
            data.append((float(x), float(y)))
    return data

# Plot the data
def plot_data(input_files, output_dir):
    plt.figure(figsize=(10, 6))
    for file_path, label in input_files:
        data = read_data(file_path)
        x, y = zip(*data)
        plt.plot(x, y, label=label)
    
    plt.xlabel(r"$\epsilon$")
    plt.ylabel(r"$\sigma$")
    plt.title(r"Constitutive Relation with Various $\alpha$")
    plt.legend()
    plt.xlim(0)
    plt.ylim(0)

    # Get the current script name and replace the extension with .png
    script_name = os.path.basename(__file__)
    output_file = os.path.join(output_dir, script_name.replace(".py", ".png"))
    plt.savefig(output_file, dpi=350)

# Define the output directory
output_directory = "plots/"

# Call the plot function with the output directory
plot_data(input_files, output_directory)