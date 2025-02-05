from .fmodule import *
from .config.load import config
#from .config import input_variables

import re
import subprocess

file = config.load_config()
load_c_config(file)

def run(cfg_file):
    folder = "/home/g/Research/fcode/build/bin/"
    #result = subprocess.run([folder + "fcode.x", cfg_file], capture_output=True, text=True)
    with open(cfg_file, 'r') as f:
        result = subprocess.run([folder + "fcode.x"], stdin=f, capture_output=True, text=True)
    return result.stdout

def forge(cfg_file, config):
    folder = "/home/g/Research/fcode/build/bin/"
    with open(cfg_file, 'w') as f:
        for section, values in config.items():
            f.write(f"[{section}]\n")
            for key, value in values.items():
                # Handle lists (e.g., k_mesh, q_mesh)
                if isinstance(value, list):
                    if all(isinstance(i, list) for i in value):  # Handle 2D lists like CELL
                        for row in value:
                            f.write("    " + " ".join(map(str, row)) + "\n")
                    else:
                        f.write(f"    {key} = " + " ".join(map(str, value)) + "\n")
                # Handle booleans (convert True/False to lowercase)
                elif isinstance(value, bool):
                    f.write(f"    {key} = {str(value).lower()}\n")
                # Handle other config types
                else:
                    f.write(f"    {key} = '{value}'\n" if isinstance(value, str) else f"    {key} = {value}\n")
            f.write("\n")  # Separate sections

def launch(config):
    cfg_file = "temp.cfg"
    forge(cfg_file, config)
    result = run(cfg_file)
    subprocess.run(["rm", cfg_file])
    return result

def grep(output, phrase):
    match = re.search(r".*?" + re.escape(phrase) + r".*", output)
    return match.group() if match else None

def extract_value(string):
    match = re.search(r"[-+]?\d*\.?\d+", string)  # Regex for both integers and floats
    return float(match.group()) 
