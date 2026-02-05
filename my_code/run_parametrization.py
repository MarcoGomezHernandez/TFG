import os
import subprocess
import glob

# Path to your C++ executable
exe_path = "./programa"

# Folder containing YAML files
yaml_folder = "yaml_configs"

# Simulation parameters
simulation_time = 20000.0  # adjust as needed
step = 0.1                # optional, adjust if needed

csv_file = os.path.abspath("PDtracesV.csv")

# Find all YAML files
yaml_files = glob.glob(os.path.join(yaml_folder, "*.yaml"))

for yaml_file in yaml_files:
    # Generate output filename by replacing .yaml with .asc
    base_name = os.path.splitext(os.path.basename(yaml_file))[0]
    output_file = os.path.join(yaml_folder, base_name + ".asc")

    # Build the command
    cmd = [exe_path, yaml_file, output_file, str(simulation_time), str(step), csv_file]

    # Run the command
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    # TODO run also python plot.py to have record of files (change eps to png in plot.py)

    # Print stdout and stderr for debugging
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)



