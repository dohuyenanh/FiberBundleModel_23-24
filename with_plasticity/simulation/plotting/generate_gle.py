import argparse
import read_parameters # pylint: disable=import-error

# Define the JSON file path
json_file = '../program/parameters.json'

# Read the required parameters from the JSON file
keys = ['E', 'e_y', 'e_b', 'de_y', 'de_b_1', 'alpha_1', 'limit']
try:
    params = read_parameters.read_parameters(json_file, *keys)
except KeyError as e:
    print(f"Error: Missing key {e} in the JSON file.")
    exit(1)

# Generate the GLE script
gle_script = f"""!Given:
E = {params['E']}           !elastic modulus of intact fibers
e_y = {params['e_y']}      !mean value of yielding thresholds
e_b = {params['e_b']}      !mean value of broken thresholds
de_y = {params['de_y']}      !radius of the uniform distribution of e_y
de_b = {params['de_b_1']}      !radius of the uniform distribution of e_b
a = {params['alpha_1']}         !ratio of the elastic modulus of yielding fibers to E
limit = {params['limit']}       !upper limit of the strain axis

!To make it easier to graph, introduce some markers at the borders between each 2 adjacent intervals:
b1 = e_y - de_y     !start yielding 
b2 = e_y + de_y     !stop intact
b3 = e_b - de_b     !start breaking
b4 = e_b + de_b     !every fiber is broken

"""

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate GLE script from JSON parameters.')
parser.add_argument('output_file', type=str, help='The output GLE file path')
args = parser.parse_args()

# Save the GLE script to the specified file
with open(args.output_file, 'w') as file:
    file.write(gle_script)

print(f"GLE script generated successfully at {args.output_file}.")