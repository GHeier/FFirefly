import ctypes
from ctypes import c_void_p, c_char_p, c_float, c_int, POINTER
import numpy as np
from pathlib import Path
import os

current_file_path = os.path.abspath(__file__)
current_file_path = current_file_path[:-47] + "build/lib/libfly.so"

lib = ctypes.CDLL(current_file_path)


class Vec(ctypes.Structure):
    _fields_ = [
        ("x", ctypes.c_float),
        ("y", ctypes.c_float),
        ("z", ctypes.c_float),
        ("w", ctypes.c_float),
        ("area", ctypes.c_float),
        ("dimension", ctypes.c_int),
        ("n", ctypes.c_int),
    ]


CALLBACKFUNC = ctypes.CFUNCTYPE(ctypes.c_float, Vec)


class Surface:
    faces: list[list[float]] = []

    def __init__(self, func=None, s_val=None):
        if func is not None and s_val is not None:
            self._callback = CALLBACKFUNC(func)  # keep reference alive
            lib.Surface_export0.argtypes = [CALLBACKFUNC, ctypes.c_float]
            lib.Surface_export0.restype = ctypes.c_void_p
            self.ptr = lib.Surface_export0(self._callback, ctypes.c_float(s_val))
        else:
            raise ValueError("Must provide func and s_val")

        if not self.ptr:
            raise RuntimeError("Failed to initialize Surface")

        # Load 'faces' field from C++
        count = ctypes.c_int()
        lib.Surface_num_faces_export0.argtypes = [ctypes.c_void_p]
        lib.Surface_num_faces_export0.restype = ctypes.c_int

        n = lib.Surface_num_faces_export0(self.ptr)
        lens = (ctypes.c_int * n)()
        total_len = 3 * len(lens)
        buf = (ctypes.c_float * total_len)()

        lib.Surface_var_faces_export0.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.Surface_var_faces_export0.restype = None
        lib.Surface_var_faces_export0(self.ptr, buf, lens, ctypes.c_int(n))

        offset = 0
        result = []
        for i in range(n):
            sub = [buf[offset + j] for j in range(lens[i])]
            result.append(sub)
            offset += lens[i]
        self.faces = result

    def __call__(self, *args):
        raise TypeError("Invalid arguments to __call__")

    def __del__(self):
        try:
            lib.destroy.argtypes = [ctypes.c_void_p]
            lib.destroy(self.ptr)
        except AttributeError:
            pass


# Begin Functions
lib.string_to_vec_export0.argtypes = [ctypes.c_char_p]
lib.string_to_vec_export0.restype = ctypes.c_void_p

def string_to_vec(arg0: str) -> any:
    result = lib.string_to_vec_export0(arg0)
    return result

lib.unpack_string_export0.argtypes = [ctypes.c_char_p]
lib.unpack_string_export0.restype = ctypes.c_void_p

def unpack_string(arg0: str) -> any:
    result = lib.unpack_string_export0(arg0)
    return result

lib.vec_to_string_export0.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_int]
lib.vec_to_string_export0.restype = ctypes.c_char_p

def vec_to_string(arg0: list[float]) -> any:
    arg0_array = np.array(arg0, dtype=np.float32)
    result = lib.vec_to_string_export0(arg0_array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), len(arg0))
    return result

lib.round_export0.argtypes = [ctypes.c_int]
lib.round_export0.restype = ctypes.c_void_p

def round(arg0: int) -> any:
    result = lib.round_export0(arg0)
    return result

lib.norm_export0.argtypes = []
lib.norm_export0.restype = ctypes.c_float

def norm() -> float:
    result = lib.norm_export0()
    return result

lib.Vec_export0.argtypes = []
lib.Vec_export0.restype = ctypes.c_void_p
lib.Vec_export1.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, ctypes.c_int]
lib.Vec_export1.restype = ctypes.c_void_p
lib.Vec_export2.argtypes = [ctypes.c_void_p]
lib.Vec_export2.restype = ctypes.c_void_p
lib.Vec_export3.argtypes = [ctypes.c_float, ctypes.c_int]
lib.Vec_export3.restype = ctypes.c_void_p

lib.Vec_x_export0.argtypes = [ctypes.c_void_p]
lib.Vec_x_export0.restype = ctypes.c_float
lib.Vec_y_export0.argtypes = [ctypes.c_void_p]
lib.Vec_y_export0.restype = ctypes.c_float
lib.Vec_z_export0.argtypes = [ctypes.c_void_p]
lib.Vec_z_export0.restype = ctypes.c_float
lib.Vec_w_export0.argtypes = [ctypes.c_void_p]
lib.Vec_w_export0.restype = ctypes.c_float
lib.Vec_area_export0.argtypes = [ctypes.c_void_p]
lib.Vec_area_export0.restype = ctypes.c_float
lib.Vec_dimension_export0.argtypes = [ctypes.c_void_p]
lib.Vec_dimension_export0.restype = ctypes.c_int
lib.Vec_n_export0.argtypes = [ctypes.c_void_p]
lib.Vec_n_export0.restype = ctypes.c_int

class Vec:
    def __init__(self, *args):
        if len(args) == 0:
            self.ptr = lib.Vec_export0()
            if not self.ptr:
                raise RuntimeError('Failed to initialize Vec')

        elif len(args) > 0 and isinstance(args[0], float):
            a = args[0]
            b = args[1] if len(args) > 1 else 0
            c = args[2] if len(args) > 2 else 0
            d = args[3] if len(args) > 3 else 0
            e = args[4] if len(args) > 4 else 0
            f = args[5] if len(args) > 5 else 3
            g = args[6] if len(args) > 6 else 1
            self.ptr = lib.Vec_export1(a, b, c, d, e, f, g)
            if not self.ptr:
                raise RuntimeError('Failed to initialize Vec')

        self.x = lib.Vec_x_export0(self.ptr)
        self.y = lib.Vec_y_export0(self.ptr)
        self.z = lib.Vec_z_export0(self.ptr)
        self.w = lib.Vec_w_export0(self.ptr)
        self.area = lib.Vec_area_export0(self.ptr)
        self.dimension = lib.Vec_dimension_export0(self.ptr)
        self.n = lib.Vec_n_export0(self.ptr)
    def __call__(self, *args):
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
