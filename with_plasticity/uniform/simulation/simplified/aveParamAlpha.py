import os

# Define the directory and output file
directory = '/Users/huyenanhdo/Library/Mobile Documents/com~apple~CloudDocs/FBM/with_plasticity/simulation/simplified/data'
output_file = 'aveParamAlpha_100000.txt'

# List to store all lines
all_lines = []

# Iterate over all files in the directory
for filename in os.listdir(directory):
    # Check if the file starts with "aveParamAlpha_"
    if filename.startswith('aveParamAlpha_100000_') and filename != 'aveParamAlpha.txt':
        file_path = os.path.join(directory, filename)
        # Open and read the content of the file
        with open(file_path, 'r') as infile:
            content = infile.readlines()
            # Add the lines to the list
            all_lines.extend(content)

# Sort the lines based on the values in the first column
all_lines.sort(key=lambda line: float(line.split()[0]))

# Open the output file in write mode and write the sorted lines
with open(output_file, 'w') as outfile:
    outfile.writelines(all_lines)