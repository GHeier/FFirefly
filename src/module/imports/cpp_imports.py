import ctypes
from ctypes import c_void_p, c_char_p, c_float, c_int, POINTER
import numpy as np
from pathlib import Path
import os

current_file_path = os.path.abspath(__file__)
current_file_path = current_file_path[:-47] + "build/lib/libfly.so"

lib = ctypes.CDLL(current_file_path)

# Begin Functions
lib.epsilon_export0.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_int]
lib.epsilon_export0.restype = ctypes.c_float

def epsilon(arg0: int, arg1: list[float]) -> float:
    arg1_array = np.array(arg1, dtype=np.float32)
    result = lib.epsilon_export0(arg0, arg1_array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(arg1))
    return result

lib.Bands_export0.argtypes = []
lib.Bands_export0.restype = ctypes.c_void_p

class Bands:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Bands_export0()
        if not self.ptr:
            raise RuntimeError('Failed to initialize Bands')

    def __call__(self, *args):
        raise TypeError('Invalid arguments to __call__')

    def __del__(self):
        try:
            destroy = lib.destroy
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass

lib.Vertex_export0.argtypes = []
lib.Vertex_export0.restype = ctypes.c_void_p

class Vertex:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Vertex_export0()
        if not self.ptr:
            raise RuntimeError('Failed to initialize Vertex')

    def __call__(self, *args):
        raise TypeError('Invalid arguments to __call__')

    def __del__(self):
        try:
            destroy = lib.destroy
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass

lib.Field_R_export0.argtypes = []
lib.Field_R_export0.restype = ctypes.c_void_p
lib.Field_R_export1.argtypes = [ctypes.c_char_p]
lib.Field_R_export1.restype = ctypes.c_void_p

lib.Field_R_operator_export0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_float]
lib.Field_R_operator_export0.restype = ctypes.c_float
lib.Field_R_operator_export1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_float]
lib.Field_R_operator_export1.restype = ctypes.c_float
lib.Field_R_operator_export2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float]
lib.Field_R_operator_export2.restype = ctypes.c_float
lib.Field_R_operator_export3.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float]
lib.Field_R_operator_export3.restype = ctypes.c_float

class Field_R:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Field_R_export0()
        else:
            self.ptr = lib.Field_R_export1(ctypes.c_char_p(filename.encode('utf-8')))
        if not self.ptr:
            raise RuntimeError('Failed to initialize Field_R')

    def __call__(self, *args):
        if len(args) == 2 and True and isinstance(args[1], (int, float)):
        arg0 = ctypes.c_float(args[1])
            return lib.Field_R_operator_export0(self.ptr, arg0)
        if len(args) == 3 and True and isinstance(args[1], int) and isinstance(args[2], (int, float)):
        arg0 = ctypes.c_int(args[1])
        arg1 = ctypes.c_float(args[2])
            return lib.Field_R_operator_export1(self.ptr, arg0, arg1)
        if len(args) == 3 and True and isinstance(args[1], (list, tuple)) and isinstance(args[2], (int, float)):
        arg0 = (ctypes.c_float * len(args[1]))(*[float(x) for x in args[1]])
        arg1 = ctypes.c_int(len(args[1]))
        arg2 = ctypes.c_float(args[2])
            return lib.Field_R_operator_export2(self.ptr, arg0, arg1, arg2)
        if len(args) == 4 and True and isinstance(args[1], int) and isinstance(args[2], (list, tuple)) and isinstance(args[3], (int, float)):
        arg0 = ctypes.c_int(args[1])
        arg1 = (ctypes.c_float * len(args[2]))(*[float(x) for x in args[2]])
        arg2 = ctypes.c_int(len(args[2]))
        arg3 = ctypes.c_float(args[3])
            return lib.Field_R_operator_export3(self.ptr, arg0, arg1, arg2, arg3)
        raise TypeError('Invalid arguments to __call__')

    def __del__(self):
        try:
            destroy = lib.destroy
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass

lib.Field_C_export0.argtypes = []
lib.Field_C_export0.restype = ctypes.c_void_p
lib.Field_C_export1.argtypes = [ctypes.c_char_p]
lib.Field_C_export1.restype = ctypes.c_void_p

lib.Field_C_operator_export0.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_float]
lib.Field_C_operator_export0.restype = ctypes.c_float
lib.Field_C_operator_export1.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_float]
lib.Field_C_operator_export1.restype = ctypes.c_float
lib.Field_C_operator_export2.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float]
lib.Field_C_operator_export2.restype = ctypes.c_float
lib.Field_C_operator_export3.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float]
lib.Field_C_operator_export3.restype = ctypes.c_float

class Field_C:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Field_C_export0()
        else:
            self.ptr = lib.Field_C_export1(ctypes.c_char_p(filename.encode('utf-8')))
        if not self.ptr:
            raise RuntimeError('Failed to initialize Field_C')

    def __call__(self, *args):
        if len(args) == 2 and True and isinstance(args[1], (int, float)):
        arg0 = ctypes.c_float(args[1])
            return lib.Field_C_operator_export0(self.ptr, arg0)
        if len(args) == 3 and True and isinstance(args[1], int) and isinstance(args[2], (int, float)):
        arg0 = ctypes.c_int(args[1])
        arg1 = ctypes.c_float(args[2])
            return lib.Field_C_operator_export1(self.ptr, arg0, arg1)
        if len(args) == 3 and True and isinstance(args[1], (list, tuple)) and isinstance(args[2], (int, float)):
        arg0 = (ctypes.c_float * len(args[1]))(*[float(x) for x in args[1]])
        arg1 = ctypes.c_int(len(args[1]))
        arg2 = ctypes.c_float(args[2])
            return lib.Field_C_operator_export2(self.ptr, arg0, arg1, arg2)
        if len(args) == 4 and True and isinstance(args[1], int) and isinstance(args[2], (list, tuple)) and isinstance(args[3], (int, float)):
        arg0 = ctypes.c_int(args[1])
        arg1 = (ctypes.c_float * len(args[2]))(*[float(x) for x in args[2]])
        arg2 = ctypes.c_int(len(args[2]))
        arg3 = ctypes.c_float(args[3])
            return lib.Field_C_operator_export3(self.ptr, arg0, arg1, arg2, arg3)
        raise TypeError('Invalid arguments to __call__')

    def __del__(self):
        try:
            destroy = lib.destroy
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass

# End Functions

# Set up the function signature
lib.load_config_export0.argtypes = [ctypes.c_char_p]
lib.load_config_export0.restype = None  # Equivalent to Cvoid


# Define the wrapper function
def load_config(path: str) -> None:
    lib.load_config_export0(path.encode("utf-8"))


# field = Field_R("sample_bands.dat")
# print(field(1, [-0.9, -0.9]))
#
# load_config("/home/g/Research/ffirefly/build/bin/input.cfg")
# print(epsilon(1, [0.1, 0.2, 0.3]))
