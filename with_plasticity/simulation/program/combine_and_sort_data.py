import os
import glob

# Define the directories to search for files
directories = ["../data/results_1", "../data/results_2"]

# Output file name
output_file = "../data/aveParamAlpha.txt"

# List to store all rows of data
all_data = []

# Search for files in the specified directories
for directory in directories:
    # Use glob to find files matching the pattern
    files = glob.glob(os.path.join(directory, "aveParamAlpha*.txt"))
    for file_path in sorted(files):  # Sort files for consistent processing
        print(f"Processing file: {file_path}")  # Optional: Show progress
        with open(file_path, "r") as infile:
            for line in infile:
                # Skip empty lines or lines that are just whitespace
                if line.strip():
                    # Split the line into columns and store as a tuple
                    columns = line.split()
                    if len(columns) >= 3:  # Ensure there are at least 3 columns
                        all_data.append((float(columns[0]), float(columns[1]), float(columns[2])))

# Sort the combined data by the first column
all_data.sort(key=lambda x: x[0])

# Write the sorted data to the output file
with open(output_file, "w") as outfile:
    for row in all_data:
        # Write each row as a space-separated line
        outfile.write(f"{row[0]} {row[1]} {row[2]}\n")

print(f"All matching files have been combined and sorted into: {output_file}")