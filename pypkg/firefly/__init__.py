from .src.config.load import config as config
from .src.module.imports.cpp_imports import *
from .plot.fly_plot import sketch

import os
import re
import subprocess

file = config.load_config()
if os.path.exists(file):
    load_config(file)
else:
    printv("No .cfg file found. Run fly.x to load in values")

current_file = os.path.abspath(__file__)
folder = current_file[:-25] + "build/bin/"


def run(cfg_file):
    # result = subprocess.run([folder + "fly.x", cfg_file], capture_output=True, text=True)
    with open(cfg_file, "r") as f:
        result = subprocess.run(
            [folder + "fly.x"], stdin=f, capture_output=True, text=True
        )
    return result.stdout, result.stderr


def forge(cfg_file, config):
    with open(cfg_file, "w") as f:
        for section, values in config.items():
            f.write(f"[{section}]\n")
            for key, value in values.items():
                # Handle lists (e.g., k_mesh, q_mesh)
                if isinstance(value, list):
                    if all(
                        isinstance(i, list) for i in value
                    ):  # Handle 2D lists like CELL
                        for row in value:
                            f.write("    " + " ".join(map(str, row)) + "\n")
                    else:
                        f.write(f"    {key} = " + " ".join(map(str, value)) + "\n")
                # Handle booleans (convert True/False to lowercase)
                elif isinstance(value, bool):
                    f.write(f"    {key} = {str(value).lower()}\n")
                # Handle other config types
                else:
                    f.write(
                        f"    {key} = '{value}'\n"
                        if isinstance(value, str)
                        else f"    {key} = {value}\n"
                    )
            f.write("\n")  # Separate sections


def launch(config):
    cfg_file = "temp.cfg"
    forge(cfg_file, config)
    output, error = run(cfg_file)
    if os.path.exists(cfg_file):
        subprocess.run(["rm", cfg_file])
    return output, error


def grep(output, phrase):
    match = re.search(r".*?" + re.escape(phrase) + r".*", output)
    return match.group() if match else None


def extract_value(string):
    match = re.search(r"[-+]?\d*\.?\d+", string)  # Regex for both integers and floats
    return float(match.group())


def unpack(filename):
    config = {}
    current_section = None
    collecting_matrix = False
    matrix_buffer = []

    def parse_value(value):
        value = value.strip().strip("'\"")
        if value.lower() in {"true", "false"}:
            return value.lower() == "true"
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        return value

    def parse_list_value(value):
        return [int(x) for x in value.split()]

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Section header
            if line.startswith("[") and line.endswith("]"):
                if collecting_matrix and current_section:
                    config[current_section]["lattice"] = matrix_buffer
                    matrix_buffer = []
                    collecting_matrix = False

                section = line[1:-1].strip()
                config[section] = {}
                current_section = section
                continue

            # Matrix-style lines (no "=")
            if "=" not in line and current_section:
                parts = line.split()
                if all(p.replace('.', '', 1).replace('-', '', 1).isdigit() for p in parts):
                    if not collecting_matrix:
                        collecting_matrix = True
                        matrix_buffer = []
                    matrix_buffer.append([float(p) for p in parts])
                    continue
                else:
                    collecting_matrix = False

            if collecting_matrix and current_section:
                config[current_section]["lattice"] = matrix_buffer
                collecting_matrix = False
                matrix_buffer = []

            # Key-value line
            if "=" in line:
                key, value = map(str.strip, line.split("=", 1))
                if key in ("k_mesh", "q_mesh"):
                    config[current_section][key] = parse_list_value(value)
                else:
                    config[current_section][key] = parse_value(value)

    # Finalize trailing matrix section
    if collecting_matrix and current_section:
        config[current_section]["lattice"] = matrix_buffer

    return config
