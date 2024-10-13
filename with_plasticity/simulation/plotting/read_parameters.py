import json

# Function to read specific parameters from a JSON file
def read_parameters(json_file, *keys):
    with open(json_file, "r", encoding="utf-8") as file:
        parameters = json.load(file)
    values = {key: parameters["values"][key] for key in keys}
    return values

# Function to read the number of specific parameters from a JSON file
def read_number_of_parameters(json_file, key):
    with open(json_file, "r", encoding="utf-8") as file:
        parameters = json.load(file)
    return parameters["number of values"][key]