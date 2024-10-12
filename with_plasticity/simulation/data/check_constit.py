def check_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    previous_value = None
    for i, line in enumerate(lines):
        columns = line.strip().split('\t')
        if len(columns) < 2:
            continue

        current_value = float(columns[1])
        if previous_value is not None and current_value < previous_value:
            print(f"Discrepancy found at line {i + 1}: {current_value} is smaller than {previous_value}")

        previous_value = current_value

# Check the file
check_file("constit.txt")