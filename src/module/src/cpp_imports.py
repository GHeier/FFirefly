import ctypes
import numpy as np
from pathlib import Path

# Define the path to the shared library
testmod = Path("/home/g/Research/fcode/build/lib/libfly.so")

# Load the shared library
lib = ctypes.CDLL(str(testmod))

# Begin Functions

# Epsilon Function
lib.epsilon_export.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_int]
lib.epsilon_export.restype = ctypes.c_float

def epsilon(n: int, k: list[float]) -> float:
    k_array = np.array(k, dtype=np.float32)
    result = lib.epsilon_export(n, k_array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(k))
    return result

# End Functions

# Set up the function signature
lib.load_config_julia.argtypes = [ctypes.c_char_p]
lib.load_config_julia.restype = None  # Equivalent to Cvoid

# Define the wrapper function
def load_config(path: str) -> None:
    lib.load_config_julia(path.encode('utf-8'))

load_config("/home/g/Research/fcode/build/bin/input.cfg")
print(epsilon(1, [0.1, 0.2, 0.3]))
