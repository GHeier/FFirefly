import ctypes
from pathlib import Path

# Define the shared library
lib_path = Path("/home/g/Research/ffirefly/build/lib/libfly.so")
lib = ctypes.CDLL(str(lib_path))

# Mapping from string to ctypes type
type_map = {
    "c_int": ctypes.c_int,
    "c_float": ctypes.c_float,
    "c_double": ctypes.c_double,
    "c_char_p": ctypes.c_char_p,
    "None": None,
    "POINTER(c_float)": ctypes.POINTER(ctypes.c_float),
    "POINTER(c_double)": ctypes.POINTER(ctypes.c_double),
}

# List of function definitions
functions = [
    {
        "name": "load_config_julia",
        "args": ["c_char_p"],
        "return": "None"
    },
    {
        "name": "epsilon_julia",
        "args": ["c_int", "POINTER(c_float)", "c_int"],
        "return": "c_float"
    }
]

# Generate wrappers
for f in functions:
    func = getattr(lib, f["name"])
    func.argtypes = [type_map[arg] for arg in f["args"]]
    func.restype = type_map[f["return"]]
    globals()[f["name"]] = func

def generate_python_function(name, arg_names):
    args_str = ", ".join(arg_names)
    call_args = ", ".join(arg_names)
    return f"""
def {name}_py({args_str}):
    return {name}({call_args})
"""

# Example usage
for f in functions:
    arg_names = [f"arg{i}" for i in range(len(f["args"]))]
    print(generate_python_function(f["name"], arg_names))
