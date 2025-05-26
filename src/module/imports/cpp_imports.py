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

lib.Field_C_export0.argtypes = []
lib.Field_C_export0.restype = ctypes.c_void_p
lib.Field_C_export1.argtypes = [ctypes.c_char_p]
lib.Field_C_export1.restype = ctypes.c_void_p
lib.Field_C_export2.restype = c_void_p

lib.Field_C_operator_export0.argtypes = [ctypes.c_void_p, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
lib.Field_C_operator_export0.restype = None
lib.Field_C_operator_export1.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
lib.Field_C_operator_export1.restype = None
lib.Field_C_operator_export2.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
lib.Field_C_operator_export2.restype = None
lib.Field_C_operator_export3.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
lib.Field_C_operator_export3.restype = None

class Field_C:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Field_C_export0()
        else:
            self.ptr = lib.Field_C_export2(c_char_p(filename.encode('utf-8')))
        if not self.ptr:
            raise RuntimeError('Failed to initialize Field_C')

    def __call__(self, *args):
        # Overload for args=1, required=1
        if len(args) >= 1 and len(args) <= 1 and isinstance(args[0], (int, float)):
            real = ctypes.c_float()
            imag = ctypes.c_float()
            arg0 = ctypes.c_float(args[0])
            lib.Field_C_operator_export0(self.ptr, arg0, ctypes.byref(real), ctypes.byref(imag))
            return complex(real.value, imag.value)
        # Overload for args=2, required=2
        if len(args) >= 2 and len(args) <= 2 and isinstance(args[0], int) and isinstance(args[1], (int, float)):
            real = ctypes.c_float()
            imag = ctypes.c_float()
            arg0 = ctypes.c_int(args[0])
            arg1 = ctypes.c_float(args[1])
            lib.Field_C_operator_export1(self.ptr, arg0, arg1, ctypes.byref(real), ctypes.byref(imag))
            return complex(real.value, imag.value)
        # Overload for args=2, required=1
        if len(args) >= 1 and len(args) <= 2:
            real = ctypes.c_float()
            imag = ctypes.c_float()
            arg0 = (ctypes.c_float * len(args[0]))(*[float(x) for x in args[0]])
            arg0_len = ctypes.c_int(len(args[0]))
            arg2 = ctypes.c_float(args[1]) if len(args) > 1 else ctypes.c_float(0.0)
            lib.Field_C_operator_export2(self.ptr, arg0, arg0_len, arg2, ctypes.byref(real), ctypes.byref(imag))
            return complex(real.value, imag.value)
        # Overload for args=3, required=2
        if len(args) >= 2 and len(args) <= 3 and isinstance(args[0], int):
            real = ctypes.c_float()
            imag = ctypes.c_float()
            arg0 = ctypes.c_int(args[0])
            arg1 = (ctypes.c_float * len(args[1]))(*[float(x) for x in args[1]])
            arg1_len = ctypes.c_int(len(args[1]))
            arg3 = ctypes.c_float(args[2]) if len(args) > 2 else ctypes.c_float(0.0)
            lib.Field_C_operator_export3(self.ptr, arg0, arg1, arg1_len, arg3, ctypes.byref(real), ctypes.byref(imag))
            return complex(real.value, imag.value)
        raise TypeError('Invalid arguments to __call__')

    def __del__(self):
        try:
            destroy = lib.destroy
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass

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
