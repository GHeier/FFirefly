import ctypes
from ctypes import c_void_p, c_char_p, c_float, c_int, POINTER
import numpy as np
from pathlib import Path

lib = ctypes.CDLL("/home/g/Research/ffirefly/build/lib/libfly.so")

# Begin Functions

# Epsilon Function
lib.epsilon_export0.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_int]
lib.epsilon_export0.restype = ctypes.c_float

def epsilon(n: int, k: list[float]) -> float:
    k_array = np.array(k, dtype=np.float32)
    result = lib.epsilon_export0(n, k_array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(k))
    return result

# End Functions

# Set up the function signature
lib.load_config_export0.argtypes = [ctypes.c_char_p]
lib.load_config_export0.restype = None  # Equivalent to Cvoid

# Define the wrapper function
def load_config(path: str) -> None:
    lib.load_config_export0(path.encode('utf-8'))

# Begin Objects

# Field_R Object
class Field_R:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Field_R_export0()
        else:
            self.ptr = lib.Field_R_export2(c_char_p(filename.encode('utf-8')))
        if not self.ptr:
            raise RuntimeError("Failed to initialize Field_R")

    def __call__(self, *args):
        # Overload for (w: float)
        if len(args) == 1 and isinstance(args[0], (int, float)):
            w = c_float(args[0])
            return lib.Field_R_operator_export0(self.ptr, w)

        # Overload for (n: int, w: float)
        if len(args) == 2 and isinstance(args[0], int) and isinstance(args[1], float):
            n = c_int(args[0])
            w = c_float(args[1])
            return lib.Field_R_operator_export1(self.ptr, n, w)

        # Overload for (k: list[float], w=0.0)
        if len(args) >= 1 and isinstance(args[0], (list, tuple)):
            k = (c_float * len(args[0]))(*[float(v) for v in args[0]])
            len_k = c_int(len(args[0]))
            w = c_float(args[1]) if len(args) == 2 else c_float(0.0)
            return lib.Field_R_operator_export2(self.ptr, k, len_k, w)

        # Overload for (n: int, k: list[float], w=0.0)
        if len(args) >= 2 and isinstance(args[0], int) and isinstance(args[1], (list, tuple)):
            n = c_int(args[0])
            k = (c_float * len(args[1]))(*[float(v) for v in args[1]])
            len_k = c_int(len(args[1]))
            w = c_float(args[2]) if len(args) == 3 else c_float(0.0)
            return lib.Field_R_operator_export3(self.ptr, n, k, len_k, w)

        raise TypeError("Invalid arguments to Field_R.__call__")

    def __del__(self):
        try:
            destroy = lib.destroy
            destroy.argtypes = [c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass  # No destructor available

# Setup return and argument types
lib.Field_R_export0.restype = c_void_p
lib.Field_R_export2.argtypes = [c_char_p]
lib.Field_R_export2.restype = c_void_p

lib.Field_R_operator_export0.argtypes = [c_void_p, c_float]
lib.Field_R_operator_export0.restype = c_float

lib.Field_R_operator_export1.argtypes = [c_void_p, c_int, c_float]
lib.Field_R_operator_export1.restype = c_float

lib.Field_R_operator_export2.argtypes = [c_void_p, POINTER(c_float), c_int, c_float]
lib.Field_R_operator_export2.restype = c_float

lib.Field_R_operator_export3.argtypes = [c_void_p, c_int, POINTER(c_float), c_int, c_float]
lib.Field_R_operator_export3.restype = c_float

#End Objects

#field = Field_R("sample_bands.dat")
#print(field(1, [-0.9, -0.9]))
#
#load_config("/home/g/Research/ffirefly/build/bin/input.cfg")
#print(epsilon(1, [0.1, 0.2, 0.3]))
