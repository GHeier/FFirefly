import ctypes
from ctypes import c_bool, c_void_p, c_char_p, c_float, c_int, POINTER, byref
import numpy as np
from pathlib import Path
import os

# Get the project root directory (FFirefly/)
current_file_path = os.path.abspath(__file__)
lib_path = current_file_path.split("FFirefly")[0] + "FFirefly/build/lib/libfly.so"

lib = ctypes.CDLL(lib_path)


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
            destroy = lib.destroy_Surface
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
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
            destroy = lib.destroy_Vec
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass

# End Functions

class CMData:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.CMData_export0()
        else:
            self.ptr = lib.CMData_export1(filename)
        if not self.ptr:
            raise RuntimeError("Failed to initialize CMData")



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

        # Load metadata from C++
        self.dimension = lib.Field_R_get_dimension(self.ptr)

        # Load mesh
        mesh_size = lib.Field_R_get_mesh_size(self.ptr)
        if mesh_size > 0:
            mesh_buf = (c_int * mesh_size)()
            lib.Field_R_get_mesh(self.ptr, mesh_buf)
            self.mesh = [mesh_buf[i] for i in range(mesh_size)]
        else:
            self.mesh = []

        # Load domain
        domain_rows = lib.Field_R_get_domain_rows(self.ptr)
        domain_cols = lib.Field_R_get_domain_cols(self.ptr)
        if domain_rows > 0 and domain_cols > 0:
            domain_buf = (c_float * (domain_rows * domain_cols))()
            lib.Field_R_get_domain(self.ptr, domain_buf)
            self.domain = [[domain_buf[i * domain_cols + j] for j in range(domain_cols)]
                          for i in range(domain_rows)]
        else:
            self.domain = []

        # Load w_points
        w_points_size = lib.Field_R_get_w_points_size(self.ptr)
        if w_points_size > 0:
            w_points_buf = (c_float * w_points_size)()
            lib.Field_R_get_w_points(self.ptr, w_points_buf)
            self.w_points = np.array([w_points_buf[i] for i in range(w_points_size)], dtype=np.float32)
        else:
            self.w_points = np.array([], dtype=np.float32)


    def __call__(self, *args):
        # Overload for (w: float)
        if len(args) == 1 and isinstance(args[0], (int, float)):
            w = c_float(args[0])
            return lib.Field_R_operator_export0(self.ptr, w)

        ## Overload for (n: int, w: float) - COMMENTED OUT, no export
        #if len(args) == 2 and isinstance(args[0], int) and isinstance(args[1], float):
        #    n = c_int(args[0])
        #    w = c_float(args[1])
        #    return lib.Field_R_operator_export1(self.ptr, n, w)

        # Overload for (k: list[float], w=0.0) or (w_points: list[float])
        if len(args) >= 1 and isinstance(args[0], (list, tuple, np.ndarray)):
            # Check if it's a list of points (list of lists)
            if len(args[0]) > 0 and isinstance(args[0][0], (list, tuple, np.ndarray)):
                # List of points
                points = args[0]
                num_points = len(points)
                if num_points == 0:
                    return np.array([], dtype=np.float32)

                point_len = len(points[0])
                points_flat = (c_float * (num_points * point_len))()
                for i, p in enumerate(points):
                    for j, val in enumerate(p):
                        points_flat[i * point_len + j] = float(val)

                w = c_float(args[1]) if len(args) == 2 else c_float(0.0)
                output = (c_float * num_points)()

                lib.Field_R_operator_export_list(self.ptr, points_flat, c_int(num_points), c_int(point_len), w, output)
                return np.array([output[i] for i in range(num_points)], dtype=np.float32)
            else:
                # Single list of numbers - could be k-point or w-points list
                # If no spatial mesh, treat as w-points list
                # Otherwise, treat as k-point
                if len(self.mesh) == 0 and len(args) == 1:
                    # List of w-points for 0D field
                    w_points = np.array(args[0], dtype=np.float32)
                    num_w = len(w_points)
                    w_array = (c_float * num_w)(*w_points)
                    output = (c_float * num_w)()
                    lib.Field_R_operator_export_w_list(self.ptr, w_array, c_int(num_w), output)
                    return np.array([output[i] for i in range(num_w)], dtype=np.float32)
                else:
                    # Single k-point
                    k = (c_float * len(args[0]))(*[float(v) for v in args[0]])
                    len_k = c_int(len(args[0]))
                    w = c_float(args[1]) if len(args) == 2 else c_float(0.0)
                    return lib.Field_R_operator_export2(self.ptr, k, len_k, w)

        ## Overload for (n: int, k: list[float], w=0.0) - COMMENTED OUT, no export
        #if len(args) >= 2 and isinstance(args[0], int) and isinstance(args[1], (list, tuple)):
        #    n = c_int(args[0])
        #    k = (c_float * len(args[1]))(*[float(v) for v in args[1]])
        #    len_k = c_int(len(args[1]))
        #    w = c_float(args[2]) if len(args) == 3 else c_float(0.0)
        #    return lib.Field_R_operator_export3(self.ptr, n, k, len_k, w)

        raise TypeError("Invalid arguments to Field_R.__call__")

    def get_data(self):
        """Get data array reshaped in (w,k) format."""
        ptr = lib.Field_R_get_data(self.ptr)
        bd = object.__new__(BaseData)
        bd.ptr = ptr
        bd.owns_ptr = False  # This is a borrowed pointer, don't destroy it
        bd._load_metadata()
        data_array = bd.get_data()

        # Reshape from (nk*nw) to (nw, nk)
        if bd.n_indices == 2:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nk, bd.nw, bd.dim_indices, bd.dim_indices)
            return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        else:
            # Scalar: (nk*nw,) -> (nw, nk)
            return data_array.reshape(bd.nk, bd.nw).T

    def __del__(self):
        try:
            destroy = lib.destroy_Field_R
            destroy.argtypes = [c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass  # No destructor available

# Setup return and argument types
lib.Field_R_export0.restype = c_void_p
lib.Field_R_export2.argtypes = [c_char_p]
lib.Field_R_export2.restype = c_void_p

# nbnd removed - not part of Field objects anymore
#lib.Field_R_nbnd_export0.argtypes = [c_void_p]
#lib.Field_R_nbnd_export0.restype = c_int

lib.Field_R_operator_export0.argtypes = [c_void_p, c_float]
lib.Field_R_operator_export0.restype = c_float

#lib.Field_R_operator_export1.argtypes = [c_void_p, c_int, c_float]
#lib.Field_R_operator_export1.restype = c_float

lib.Field_R_operator_export2.argtypes = [c_void_p, POINTER(c_float), c_int, c_float]
lib.Field_R_operator_export2.restype = c_float

#lib.Field_R_operator_export3.argtypes = [c_void_p, c_int, POINTER(c_float), c_int, c_float]
#lib.Field_R_operator_export3.restype = c_float

lib.Field_R_operator_export_list.argtypes = [c_void_p, POINTER(c_float), c_int, c_int, c_float, POINTER(c_float)]
lib.Field_R_operator_export_list.restype = None

lib.Field_R_operator_export_w_list.argtypes = [c_void_p, POINTER(c_float), c_int, POINTER(c_float)]
lib.Field_R_operator_export_w_list.restype = None

# Field_R metadata functions
lib.Field_R_get_mesh_size.argtypes = [c_void_p]
lib.Field_R_get_mesh_size.restype = c_int
lib.Field_R_get_mesh.argtypes = [c_void_p, POINTER(c_int)]
lib.Field_R_get_mesh.restype = None
lib.Field_R_get_domain_rows.argtypes = [c_void_p]
lib.Field_R_get_domain_rows.restype = c_int
lib.Field_R_get_domain_cols.argtypes = [c_void_p]
lib.Field_R_get_domain_cols.restype = c_int
lib.Field_R_get_domain.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_R_get_domain.restype = None
lib.Field_R_get_w_points_size.argtypes = [c_void_p]
lib.Field_R_get_w_points_size.restype = c_int
lib.Field_R_get_w_points.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_R_get_w_points.restype = None
lib.Field_R_get_dimension.argtypes = [c_void_p]
lib.Field_R_get_dimension.restype = c_int

#End Objects

lib.Field_C_export0.argtypes = []
lib.Field_C_export0.restype = ctypes.c_void_p
#lib.Field_C_export1.argtypes = [ctypes.c_char_p]
#lib.Field_C_export1.restype = ctypes.c_void_p
lib.Field_C_export2.argtypes = [ctypes.c_char_p]
lib.Field_C_export2.restype = c_void_p

# nbnd removed - not part of Field objects anymore
#lib.Field_C_nbnd_export0.argtypes = [c_void_p]
#lib.Field_C_nbnd_export0.restype = c_int

lib.Field_C_operator_export0.argtypes = [ctypes.c_void_p, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
lib.Field_C_operator_export0.restype = None
#lib.Field_C_operator_export1.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
#lib.Field_C_operator_export1.restype = None
lib.Field_C_operator_export2.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
lib.Field_C_operator_export2.restype = None
#lib.Field_C_operator_export3.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float)]
#lib.Field_C_operator_export3.restype = None

lib.Field_C_operator_export_list.argtypes = [c_void_p, POINTER(c_float), c_int, c_int, c_float, POINTER(c_float), POINTER(c_float)]
lib.Field_C_operator_export_list.restype = None

lib.Field_C_operator_export_w_list.argtypes = [c_void_p, POINTER(c_float), c_int, POINTER(c_float), POINTER(c_float)]
lib.Field_C_operator_export_w_list.restype = None

# Field_C metadata functions
lib.Field_C_get_mesh_size.argtypes = [c_void_p]
lib.Field_C_get_mesh_size.restype = c_int
lib.Field_C_get_mesh.argtypes = [c_void_p, POINTER(c_int)]
lib.Field_C_get_mesh.restype = None
lib.Field_C_get_domain_rows.argtypes = [c_void_p]
lib.Field_C_get_domain_rows.restype = c_int
lib.Field_C_get_domain_cols.argtypes = [c_void_p]
lib.Field_C_get_domain_cols.restype = c_int
lib.Field_C_get_domain.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_C_get_domain.restype = None
lib.Field_C_get_w_points_size.argtypes = [c_void_p]
lib.Field_C_get_w_points_size.restype = c_int
lib.Field_C_get_w_points.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_C_get_w_points.restype = None
lib.Field_C_get_dimension.argtypes = [c_void_p]
lib.Field_C_get_dimension.restype = c_int

class Field_C:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Field_C_export0()
        else:
            self.ptr = lib.Field_C_export2(c_char_p(filename.encode('utf-8')))
        if not self.ptr:
            raise RuntimeError('Failed to initialize Field_C')

        # Load metadata from C++
        self.dimension = lib.Field_C_get_dimension(self.ptr)

        # Load mesh
        mesh_size = lib.Field_C_get_mesh_size(self.ptr)
        if mesh_size > 0:
            mesh_buf = (c_int * mesh_size)()
            lib.Field_C_get_mesh(self.ptr, mesh_buf)
            self.mesh = [mesh_buf[i] for i in range(mesh_size)]
        else:
            self.mesh = []

        # Load domain
        domain_rows = lib.Field_C_get_domain_rows(self.ptr)
        domain_cols = lib.Field_C_get_domain_cols(self.ptr)
        if domain_rows > 0 and domain_cols > 0:
            domain_buf = (c_float * (domain_rows * domain_cols))()
            lib.Field_C_get_domain(self.ptr, domain_buf)
            self.domain = [[domain_buf[i * domain_cols + j] for j in range(domain_cols)]
                          for i in range(domain_rows)]
        else:
            self.domain = []

        # Load w_points
        w_points_size = lib.Field_C_get_w_points_size(self.ptr)
        if w_points_size > 0:
            w_points_buf = (c_float * w_points_size)()
            lib.Field_C_get_w_points(self.ptr, w_points_buf)
            self.w_points = np.array([w_points_buf[i] for i in range(w_points_size)], dtype=np.float32)
        else:
            self.w_points = np.array([], dtype=np.float32)

    def __call__(self, *args):
        # Overload for (w: float)
        if len(args) == 1 and isinstance(args[0], (int, float, np.float32, np.float64)):
            real = ctypes.c_float()
            imag = ctypes.c_float()
            arg0 = ctypes.c_float(args[0])
            lib.Field_C_operator_export0(self.ptr, arg0, ctypes.byref(real), ctypes.byref(imag))
            return complex(real.value, imag.value)

        ## Overload for (n: int, w: float) - COMMENTED OUT, no export
        #if len(args) == 2 and isinstance(args[0], int) and isinstance(args[1], (int, float)):
        #    real = ctypes.c_float()
        #    imag = ctypes.c_float()
        #    arg0 = ctypes.c_int(args[0])
        #    arg1 = ctypes.c_float(args[1])
        #    lib.Field_C_operator_export1(self.ptr, arg0, arg1, ctypes.byref(real), ctypes.byref(imag))
        #    return complex(real.value, imag.value)

        # Overload for (k: list[float], w=0.0) or (w_points: list[float])
        if len(args) >= 1 and isinstance(args[0], (list, tuple, np.ndarray)):
            # Check if it's a list of points (list of lists)
            if len(args[0]) > 0 and isinstance(args[0][0], (list, tuple, np.ndarray)):
                # List of points
                points = args[0]
                num_points = len(points)
                if num_points == 0:
                    return np.array([], dtype=np.complex64)

                point_len = len(points[0])
                points_flat = (c_float * (num_points * point_len))()
                for i, p in enumerate(points):
                    for j, val in enumerate(p):
                        points_flat[i * point_len + j] = float(val)

                w = c_float(args[1]) if len(args) > 1 else c_float(0.0)
                real_output = (c_float * num_points)()
                imag_output = (c_float * num_points)()

                lib.Field_C_operator_export_list(self.ptr, points_flat, c_int(num_points), c_int(point_len), w, real_output, imag_output)
                return np.array([complex(real_output[i], imag_output[i]) for i in range(num_points)], dtype=np.complex64)
            else:
                # Single list of numbers - could be k-point or w-points list
                # If no spatial mesh, treat as w-points list
                # Otherwise, treat as k-point
                if len(self.mesh) == 0 and len(args) == 1:
                    # List of w-points for 0D field
                    w_points = np.array(args[0], dtype=np.float32)
                    num_w = len(w_points)
                    w_array = (c_float * num_w)(*w_points)
                    real_output = (c_float * num_w)()
                    imag_output = (c_float * num_w)()
                    lib.Field_C_operator_export_w_list(self.ptr, w_array, c_int(num_w), real_output, imag_output)
                    return np.array([complex(real_output[i], imag_output[i]) for i in range(num_w)], dtype=np.complex64)
                else:
                    # Single k-point
                    real = ctypes.c_float()
                    imag = ctypes.c_float()
                    arg0 = (ctypes.c_float * len(args[0]))(*[float(x) for x in args[0]])
                    arg0_len = ctypes.c_int(len(args[0]))
                    arg2 = ctypes.c_float(args[1]) if len(args) > 1 else ctypes.c_float(0.0)
                    lib.Field_C_operator_export2(self.ptr, arg0, arg0_len, arg2, ctypes.byref(real), ctypes.byref(imag))
                    return complex(real.value, imag.value)
        ## Overload for args=2-3 (n: int, k: list[float], w=0.0) - COMMENTED OUT, no export
        #if len(args) >= 2 and len(args) <= 3 and isinstance(args[0], int):
        #    real = ctypes.c_float()
        #    imag = ctypes.c_float()
        #    arg0 = ctypes.c_int(args[0])
        #    arg1 = (ctypes.c_float * len(args[1]))(*[float(x) for x in args[1]])
        #    arg1_len = ctypes.c_int(len(args[1]))
        #    arg3 = ctypes.c_float(args[2]) if len(args) > 2 else ctypes.c_float(0.0)
        #    lib.Field_C_operator_export3(self.ptr, arg0, arg1, arg1_len, arg3, ctypes.byref(real), ctypes.byref(imag))
        #    return complex(real.value, imag.value)
        raise TypeError('Invalid arguments to __call__')

    def get_data(self):
        """Get data array reshaped in (w,k) format."""
        ptr = lib.Field_C_get_data(self.ptr)
        bd = object.__new__(BaseData)
        bd.ptr = ptr
        bd.owns_ptr = False  # This is a borrowed pointer, don't destroy it
        bd._load_metadata()
        data_array = bd.get_data()

        # Reshape from (nk*nw) to (nw, nk)
        if bd.n_indices == 2:
            print("reshaping matrix data")
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nk, bd.nw, bd.dim_indices, bd.dim_indices)
            return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        else:
            print("reshaping scalar data")
            # Scalar: (nk*nw,) -> (nw, nk)
            return data_array.reshape(bd.nw, bd.nk)
            #return data_array.reshape(bd.nk, bd.nw)

    def __del__(self):
        try:
            if hasattr(self, 'ptr') and self.ptr:
                destroy = lib.destroy_Field_C
                destroy.argtypes = [ctypes.c_void_p]
                destroy(self.ptr)
                self.ptr = None
        except (AttributeError, OSError):
            pass

# Hamiltonian class
lib.Hamiltonian_export0.restype = c_void_p
lib.Hamiltonian_operator_export0.argtypes = [c_void_p, POINTER(c_float), c_int, c_float, POINTER(c_float), POINTER(c_float), POINTER(c_int)]
lib.Hamiltonian_operator_export0.restype = None

class Hamiltonian:
    def __init__(self):
        self.ptr = lib.Hamiltonian_export0()
        if not self.ptr:
            raise RuntimeError('Failed to initialize Hamiltonian')

    def __call__(self, k, w=0.0):
        """
        Evaluate Hamiltonian H(k, w) and return as numpy matrix.

        Args:
            k: momentum point (list or array)
            w: frequency (default 0.0)

        Returns:
            numpy array of shape (n, n) with complex values
        """
        k_array = (c_float * len(k))(*[float(x) for x in k])
        k_len = c_int(len(k))
        w_val = c_float(w)

        # Matrix size will be determined by the C++ function
        matrix_size = c_int(0)

        # Allocate space for a maximum size matrix (e.g., 100x100)
        max_size = 100
        real_result = (c_float * (max_size * max_size))()
        imag_result = (c_float * (max_size * max_size))()

        lib.Hamiltonian_operator_export0(
            self.ptr, k_array, k_len, w_val,
            real_result, imag_result, ctypes.byref(matrix_size)
        )

        n = matrix_size.value
        if n == 0:
            return np.array([[]], dtype=np.complex64)

        # Reshape flattened arrays to matrices
        real_matrix = np.array([real_result[i] for i in range(n*n)]).reshape(n, n)
        imag_matrix = np.array([imag_result[i] for i in range(n*n)]).reshape(n, n)

        return real_matrix + 1j * imag_matrix

    def __del__(self):
        try:
            destroy = lib.destroy_Hamiltonian
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            print("failed to clear memory")

# Field_RM class (Real Matrix field)
lib.Field_RM_export0.restype = c_void_p
lib.Field_RM_export2.argtypes = [c_char_p]
lib.Field_RM_export2.restype = c_void_p
lib.Field_RM_operator_export0.argtypes = [c_void_p, POINTER(c_float), c_int, c_float, POINTER(c_float), POINTER(c_int)]
lib.Field_RM_operator_export0.restype = None
lib.Field_RM_operator_export_list.argtypes = [c_void_p, POINTER(c_float), c_int, c_int, c_float, POINTER(c_float), POINTER(c_int)]
lib.Field_RM_operator_export_list.restype = None

# Field_RM metadata functions
lib.Field_RM_get_mesh_size.argtypes = [c_void_p]
lib.Field_RM_get_mesh_size.restype = c_int
lib.Field_RM_get_mesh.argtypes = [c_void_p, POINTER(c_int)]
lib.Field_RM_get_mesh.restype = None
lib.Field_RM_get_domain_rows.argtypes = [c_void_p]
lib.Field_RM_get_domain_rows.restype = c_int
lib.Field_RM_get_domain_cols.argtypes = [c_void_p]
lib.Field_RM_get_domain_cols.restype = c_int
lib.Field_RM_get_domain.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_RM_get_domain.restype = None
lib.Field_RM_get_w_points_size.argtypes = [c_void_p]
lib.Field_RM_get_w_points_size.restype = c_int
lib.Field_RM_get_w_points.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_RM_get_w_points.restype = None
lib.Field_RM_get_dimension.argtypes = [c_void_p]
lib.Field_RM_get_dimension.restype = c_int

class Field_RM:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Field_RM_export0()
        else:
            self.ptr = lib.Field_RM_export2(c_char_p(filename.encode('utf-8')))
        if not self.ptr:
            raise RuntimeError('Failed to initialize Field_RM')

        # Load metadata from C++
        self.dimension = lib.Field_RM_get_dimension(self.ptr)

        # Load mesh
        mesh_size = lib.Field_RM_get_mesh_size(self.ptr)
        if mesh_size > 0:
            mesh_buf = (c_int * mesh_size)()
            lib.Field_RM_get_mesh(self.ptr, mesh_buf)
            self.mesh = [mesh_buf[i] for i in range(mesh_size)]
        else:
            self.mesh = []

        # Load domain
        domain_rows = lib.Field_RM_get_domain_rows(self.ptr)
        domain_cols = lib.Field_RM_get_domain_cols(self.ptr)
        if domain_rows > 0 and domain_cols > 0:
            domain_buf = (c_float * (domain_rows * domain_cols))()
            lib.Field_RM_get_domain(self.ptr, domain_buf)
            self.domain = [[domain_buf[i * domain_cols + j] for j in range(domain_cols)]
                          for i in range(domain_rows)]
        else:
            self.domain = []

        # Load w_points
        w_points_size = lib.Field_RM_get_w_points_size(self.ptr)
        if w_points_size > 0:
            w_points_buf = (c_float * w_points_size)()
            lib.Field_RM_get_w_points(self.ptr, w_points_buf)
            self.w_points = np.array([w_points_buf[i] for i in range(w_points_size)], dtype=np.float32)
        else:
            self.w_points = np.array([], dtype=np.float32)

    def __call__(self, k, w=0.0):
        """
        Evaluate Field_RM at point(s) k with frequency w and return as numpy matrix or list of matrices.

        Args:
            k: momentum point (list or array) or list of momentum points
            w: frequency (default 0.0)

        Returns:
            numpy array of shape (n, n) with real values, or list of such arrays
        """
        # Check if k is a list of points (list of lists)
        if len(k) > 0 and isinstance(k[0], (list, tuple, np.ndarray)):
            # List of points
            points = k
            num_points = len(points)
            if num_points == 0:
                return []

            point_len = len(points[0])
            points_flat = (c_float * (num_points * point_len))()
            for i, p in enumerate(points):
                for j, val in enumerate(p):
                    points_flat[i * point_len + j] = float(val)

            w_val = c_float(w)
            matrix_size = c_int(0)

            # Allocate space for multiple matrices
            max_size = 100
            output = (c_float * (num_points * max_size * max_size))()

            lib.Field_RM_operator_export_list(
                self.ptr, points_flat, c_int(num_points), c_int(point_len), w_val,
                output, ctypes.byref(matrix_size)
            )

            n = matrix_size.value
            if n == 0:
                return [np.array([[]], dtype=np.float32) for _ in range(num_points)]

            # Reshape to list of matrices
            matrices = []
            for i in range(num_points):
                start_idx = i * n * n
                end_idx = (i + 1) * n * n
                matrix = np.array([output[j] for j in range(start_idx, end_idx)]).reshape(n, n)
                matrices.append(matrix)

            return matrices
        else:
            # Single point
            k_array = (c_float * len(k))(*[float(x) for x in k])
            k_len = c_int(len(k))
            w_val = c_float(w)

            # Matrix size will be determined by the C++ function
            matrix_size = c_int(0)

            # Allocate space for a maximum size matrix (e.g., 100x100)
            max_size = 100
            result = (c_float * (max_size * max_size))()

            lib.Field_RM_operator_export0(
                self.ptr, k_array, k_len, w_val,
                result, ctypes.byref(matrix_size)
            )

            n = matrix_size.value
            if n == 0:
                return np.array([[]], dtype=np.float32)

            # Reshape flattened array to matrix
            matrix = np.array([result[i] for i in range(n*n)]).reshape(n, n)

            return matrix

    def get_data(self):
        """Get data array reshaped in (w,k) format."""
        ptr = lib.Field_RM_get_data(self.ptr)
        bd = object.__new__(BaseData)
        bd.ptr = ptr
        bd.owns_ptr = False  # This is a borrowed pointer, don't destroy it
        bd._load_metadata()
        data_array = bd.get_data()

        # Reshape from (nk*nw) to (nw, nk)
        if bd.n_indices == 2:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nk, bd.nw, bd.dim_indices, bd.dim_indices)
            return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        else:
            # Scalar: (nk*nw,) -> (nw, nk)
            return data_array.reshape(bd.nk, bd.nw).T

    def __del__(self):
        try:
            destroy = lib.destroy_Field_RM
            destroy.argtypes = [c_void_p]
            destroy(self.ptr)
        except AttributeError:
            print("failed to clear memory")

# Field_CM class (Complex Matrix field)
lib.Field_CM_export0.restype = c_void_p
lib.Field_CM_export2.argtypes = [c_char_p]
lib.Field_CM_export2.restype = c_void_p
lib.Field_CM_operator_export0.argtypes = [c_void_p, POINTER(c_float), c_int, c_float, POINTER(c_float), POINTER(c_float), POINTER(c_int)]
lib.Field_CM_operator_export0.restype = None
lib.Field_CM_operator_export_list.argtypes = [c_void_p, POINTER(c_float), c_int, c_int, c_float, POINTER(c_float), POINTER(c_float), POINTER(c_int)]
lib.Field_CM_operator_export_list.restype = None

# Field_CM metadata functions
lib.Field_CM_get_mesh_size.argtypes = [c_void_p]
lib.Field_CM_get_mesh_size.restype = c_int
lib.Field_CM_get_mesh.argtypes = [c_void_p, POINTER(c_int)]
lib.Field_CM_get_mesh.restype = None
lib.Field_CM_get_domain_rows.argtypes = [c_void_p]
lib.Field_CM_get_domain_rows.restype = c_int
lib.Field_CM_get_domain_cols.argtypes = [c_void_p]
lib.Field_CM_get_domain_cols.restype = c_int
lib.Field_CM_get_domain.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_CM_get_domain.restype = None
lib.Field_CM_get_w_points_size.argtypes = [c_void_p]
lib.Field_CM_get_w_points_size.restype = c_int
lib.Field_CM_get_w_points.argtypes = [c_void_p, POINTER(c_float)]
lib.Field_CM_get_w_points.restype = None
lib.Field_CM_get_dimension.argtypes = [c_void_p]
lib.Field_CM_get_dimension.restype = c_int

class Field_CM:
    def __init__(self, filename=None):
        if filename is None:
            self.ptr = lib.Field_CM_export0()
        else:
            self.ptr = lib.Field_CM_export2(c_char_p(filename.encode('utf-8')))
        if not self.ptr:
            raise RuntimeError('Failed to initialize Field_CM')

        # Load metadata from C++
        self.dimension = lib.Field_CM_get_dimension(self.ptr)

        # Load mesh
        mesh_size = lib.Field_CM_get_mesh_size(self.ptr)
        if mesh_size > 0:
            mesh_buf = (c_int * mesh_size)()
            lib.Field_CM_get_mesh(self.ptr, mesh_buf)
            self.mesh = [mesh_buf[i] for i in range(mesh_size)]
        else:
            self.mesh = []

        # Load domain
        domain_rows = lib.Field_CM_get_domain_rows(self.ptr)
        domain_cols = lib.Field_CM_get_domain_cols(self.ptr)
        if domain_rows > 0 and domain_cols > 0:
            domain_buf = (c_float * (domain_rows * domain_cols))()
            lib.Field_CM_get_domain(self.ptr, domain_buf)
            self.domain = [[domain_buf[i * domain_cols + j] for j in range(domain_cols)]
                          for i in range(domain_rows)]
        else:
            self.domain = []

        # Load w_points
        w_points_size = lib.Field_CM_get_w_points_size(self.ptr)
        if w_points_size > 0:
            w_points_buf = (c_float * w_points_size)()
            lib.Field_CM_get_w_points(self.ptr, w_points_buf)
            self.w_points = np.array([w_points_buf[i] for i in range(w_points_size)], dtype=np.float32)
        else:
            self.w_points = np.array([], dtype=np.float32)

    def __call__(self, k, w=0.0):
        """
        Evaluate Field_CM at point(s) k with frequency w and return as numpy matrix or list of matrices.

        Args:
            k: momentum point (list or array) or list of momentum points
            w: frequency (default 0.0)

        Returns:
            numpy array of shape (n, n) with complex values, or list of such arrays
        """
        # Check if k is a list of points (list of lists)
        if len(k) > 0 and isinstance(k[0], (list, tuple, np.ndarray)):
            # List of points
            points = k
            num_points = len(points)
            if num_points == 0:
                return []

            point_len = len(points[0])
            points_flat = (c_float * (num_points * point_len))()
            for i, p in enumerate(points):
                for j, val in enumerate(p):
                    points_flat[i * point_len + j] = float(val)

            w_val = c_float(w)
            matrix_size = c_int(0)

            # Allocate space for multiple matrices
            max_size = 100
            real_output = (c_float * (num_points * max_size * max_size))()
            imag_output = (c_float * (num_points * max_size * max_size))()

            lib.Field_CM_operator_export_list(
                self.ptr, points_flat, c_int(num_points), c_int(point_len), w_val,
                real_output, imag_output, ctypes.byref(matrix_size)
            )

            n = matrix_size.value
            if n == 0:
                return [np.array([[]], dtype=np.complex64) for _ in range(num_points)]

            # Reshape to list of matrices
            matrices = []
            for i in range(num_points):
                start_idx = i * n * n
                end_idx = (i + 1) * n * n
                real_mat = np.array([real_output[j] for j in range(start_idx, end_idx)]).reshape(n, n)
                imag_mat = np.array([imag_output[j] for j in range(start_idx, end_idx)]).reshape(n, n)
                matrices.append(real_mat + 1j * imag_mat)

            return matrices
        else:
            # Single point
            k_array = (c_float * len(k))(*[float(x) for x in k])
            k_len = c_int(len(k))
            w_val = c_float(w)

            # Matrix size will be determined by the C++ function
            matrix_size = c_int(0)

            # Allocate space for a maximum size matrix (e.g., 100x100)
            max_size = 100
            real_result = (c_float * (max_size * max_size))()
            imag_result = (c_float * (max_size * max_size))()

            lib.Field_CM_operator_export0(
                self.ptr, k_array, k_len, w_val,
                real_result, imag_result, ctypes.byref(matrix_size)
            )

        n = matrix_size.value
        if n == 0:
            return np.array([[]], dtype=np.complex64)

        # Reshape flattened arrays to matrices
        real_matrix = np.array([real_result[i] for i in range(n*n)]).reshape(n, n)
        imag_matrix = np.array([imag_result[i] for i in range(n*n)]).reshape(n, n)

        return real_matrix + 1j * imag_matrix

    def get_data(self):
        """Get data array reshaped in (w,k) format."""
        ptr = lib.Field_CM_get_data(self.ptr)
        bd = object.__new__(BaseData)
        bd.ptr = ptr
        bd.owns_ptr = False  # This is a borrowed pointer, don't destroy it
        bd._load_metadata()
        data_array = bd.get_data()

        # Reshape from (nk*nw) to (nw, nk)
        if bd.n_indices == 4:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            dim = bd.dim_indices
            data_reshaped = data_array.reshape(bd.nw, bd.nk, dim, dim, dim, dim)
            return data_reshaped
            #return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        if bd.n_indices == 2:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nw, bd.nk, bd.dim_indices, bd.dim_indices)
            return data_reshaped
            #return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        else:
            # Scalar: (nk*nw,) -> (nw, nk)
            return data_array.reshape(bd.nw, bd.nk)

    def __del__(self):
        try:
            destroy = lib.destroy_Field_CM
            destroy.argtypes = [c_void_p]
            destroy(self.ptr)
        except AttributeError:
            print("failed to clear memory")

# Bands object
lib.Bands_export0.restype = c_void_p

lib.Bands_operator_export0.restype = c_float
lib.Bands_operator_export0.argtypes = [c_void_p, c_int, POINTER(c_float), c_int]

lib.Bands_operator_export1.restype = c_float
lib.Bands_operator_export1.argtypes = [c_void_p, POINTER(c_float), c_int]

lib.Bands_operator_export0_numpy.restype = None
lib.Bands_operator_export0_numpy.argtypes = [
    ctypes.c_void_p,            # Bands* obj
    ctypes.c_int,               # int n
    ctypes.POINTER(ctypes.c_float),  # const float* points
    ctypes.c_int,               # int num_points
    ctypes.c_int,               # int len
    ctypes.POINTER(ctypes.c_float)   # float* output
]

lib.Bands_operator_export1_numpy.restype = None
lib.Bands_operator_export1_numpy.argtypes = [
    ctypes.c_void_p,            # Bands* obj
    ctypes.POINTER(ctypes.c_float),  # const float* points
    ctypes.c_int,               # int num_points
    ctypes.c_int,               # int len
    ctypes.POINTER(ctypes.c_float)   # float* output
]

class Bands:
    def __init__(self):
        self.ptr = lib.Bands_export0()
        if not self.ptr:
            raise RuntimeError('Failed to initialize Bands')

    def __call__(self, *args):
        # Overload for (n: int, k: list/array) - single point
        if len(args) == 2 and isinstance(args[0], int) and isinstance(args[1], (list, tuple)):
            n = ctypes.c_int(args[0])
            k = (ctypes.c_float * len(args[1]))(*[float(x) for x in args[1]])
            klen = len(args[1])
            return lib.Bands_operator_export0(self.ptr, n, k, klen)

        # Overload for (n: int, points: np.ndarray) - multiple points
        if len(args) == 2 and isinstance(args[0], int) and isinstance(args[1], np.ndarray):
            kpts = args[1]
            if kpts.ndim != 2:
                raise ValueError("points must be a 2D numpy array")
            if kpts.dtype != np.float32:
                kpts = kpts.astype(np.float32)

            num_points, klen = kpts.shape

            # Allocate output array
            output = np.empty(num_points, dtype=np.float32)

            # Convert input and output to ctypes pointers
            points_ctypes = kpts.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            output_ctypes = output.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

            lib.Bands_operator_export0_numpy(
                self.ptr,
                ctypes.c_int(args[0]),
                points_ctypes,
                ctypes.c_int(num_points),
                ctypes.c_int(klen),
                output_ctypes
            )

            return output

        # Overload for (k: list/array) - single point, no band index
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            k = (ctypes.c_float * len(args[0]))(*[float(x) for x in args[0]])
            klen = len(args[0])
            return lib.Bands_operator_export1(self.ptr, k, klen)

        # Overload for (points: np.ndarray) - multiple points, no band index
        if len(args) == 1 and isinstance(args[0], np.ndarray):
            kpts = args[0]
            if kpts.ndim != 2:
                raise ValueError("points must be a 2D numpy array")
            if kpts.dtype != np.float32:
                kpts = kpts.astype(np.float32)

            num_points, klen = kpts.shape

            # Allocate output array
            output = np.empty(num_points, dtype=np.float32)

            # Convert input and output to ctypes pointers
            points_ctypes = kpts.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            output_ctypes = output.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

            lib.Bands_operator_export1_numpy(
                self.ptr,
                points_ctypes,
                ctypes.c_int(num_points),
                ctypes.c_int(klen),
                output_ctypes
            )

            return output

        raise TypeError(f"Invalid arguments to Bands.__call__: {args}")

    def __del__(self):
        try:
            destroy = lib.destroy_Bands
            destroy.argtypes = [ctypes.c_void_p]
            destroy(self.ptr)
        except AttributeError:
            pass

ctypes.POINTER(ctypes.c_float)
lib.data_save_export0.argtypes = [c_char_p, POINTER(c_float), POINTER(c_float), c_int, c_bool, c_bool, c_bool, c_bool]
lib.data_save_export0.restype = None
def data_save(filename, points, values, dimension, with_w, with_n, is_complex, is_vector):
    # Ensure inputs are numpy arrays of type float32
    points = np.asarray(points, dtype=np.float32)
    values = np.asarray(values, dtype=np.float32)

    # Calculate number of points
    if points.ndim == 2:
        num_points = points.shape[0]
    elif points.ndim == 1:
        if dimension == 0:
            raise ValueError("dimension cannot be zero when points is 1D")
        num_points = len(points) // dimension
    else:
        raise ValueError("points must be 1D or 2D array")

    # Call the C++ function
    lib.data_save_export0(
        filename.encode('utf-8'),
        points.ctypes.data_as(POINTER(c_float)),
        values.ctypes.data_as(POINTER(c_float)),
        c_int(num_points),
        c_int(dimension),
        c_bool(with_w),
        c_bool(with_n),
        c_bool(is_complex),
        c_bool(is_vector)
    )

# Define interleave_complex (similar to Julia's interleave)
def interleave_complex(values: np.ndarray) -> np.ndarray:
    values = values.astype(np.complex64)
    # For 3D arrays (k-space + frequency), use order='F' for w-k ordering
    # For 2D arrays, use default C order (row-major)
    order = 'F' if values.ndim == 3 else 'C'
    real = np.real(values).ravel(order=order)
    imag = np.imag(values).ravel(order=order)
    return np.column_stack((real, imag)).astype(np.float32).ravel()

def save_field(filename: str, values, domain, mesh, w_points=None):
    print("1")
    if w_points is None:
        w_points = []

    print("2")
    with_w = len(w_points) > 0
    found_mesh = values.shape
    nbnd = 1

    print("3")
    if with_w and found_mesh[0] != len(w_points) or (not with_w and found_mesh[0] != mesh[0]):
        nbnd = found_mesh[0]

    print("4")
    domain = np.reshape(domain, -1).astype(np.float32)
    mesh_arr = np.array(mesh, dtype=np.int32)
    len_mesh = len(mesh_arr)
    print("5")

    is_complex = np.iscomplexobj(values)
    is_vector = False  # Vector support not implemented
    with_n = nbnd > 1
    print("6")

    if not is_complex:
        values_c = np.reshape(values, -1).astype(np.float32)
    else:
        values_c = interleave_complex(values)
    print("7")

    w_points = np.array(w_points, dtype=np.float32) if with_w else np.array([], dtype=np.float32)
    print("8")

    # Set up argument types
    lib.field_save_export0.argtypes = [
        c_char_p,
        POINTER(c_float),
        POINTER(c_int),
        c_int,
        c_int,
        POINTER(c_float),
        c_int,
        c_bool,
        c_bool,
        c_bool,
        c_bool,
        POINTER(c_float),
    ]
    lib.field_save_export0.restype = None
    print("9")

    # Call C function
    print(values)
    print(filename)
    print(domain)
    print(mesh_arr)
    print(len_mesh)
    print(nbnd)
    print(w_points)
    print(len(w_points))
    print(is_complex)
    print(is_vector)
    print(with_w)
    print(with_n)
    lib.field_save_export0(
        filename.encode("utf-8"),
        domain.ctypes.data_as(POINTER(c_float)),
        mesh_arr.ctypes.data_as(POINTER(c_int)),
        len_mesh,
        nbnd,
        w_points.ctypes.data_as(POINTER(c_float)),
        len(w_points),
        is_complex,
        is_vector,
        with_w,
        with_n,
        values_c.ctypes.data_as(POINTER(c_float))
    )
    print("10")

# Set up the function signature
lib.load_config_export0.argtypes = [ctypes.c_char_p]
lib.load_config_export0.restype = None  # Equivalent to Cvoid

# Define the wrapper function
def load_config(path: str) -> None:
    lib.load_config_export0(path.encode("utf-8"))

# Save data functions
def save_data(filename: str, data: np.ndarray, mesh=None, domain=None,
              w_points=None, n_indices=1, dim_indices=1):
    """Save data to HDF5 file with automatic dispatch based on data type.

    Args:
        filename: Output filename
        data: nD array (complex or real)
        mesh: Mesh dimensions (default: empty array)
        domain: Domain vectors (default: zeros)
        w_points: Frequency points (default: empty array)
        n_indices: Number of band indices (default: 1)
        dim_indices: Dimension of indices (1=scalar, >1=matrix) (default: 1)
    """
    # Set defaults
    if mesh is None:
        mesh = np.array([], dtype=np.int32)
    if domain is None:
        domain = np.zeros((0, 0), dtype=np.float32)
    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    is_complex = np.iscomplexobj(data)

    # Determine if data is matrix or scalar (ignoring vector for now)
    if dim_indices > 1:
        # Matrix data
        # Calculate number of matrices based on data shape
        data_shape = data.shape
        num_matrices = np.prod(data_shape[:-2]) if len(data_shape) > 2 else 1
        mat_dim = dim_indices
        save_data_matrix(filename, data, num_matrices, mat_dim, is_complex,
                        mesh, domain, w_points)
    else:
        # Scalar data
        save_data_scalar(filename, data, is_complex, mesh, domain, w_points)

lib.save_data_scalar_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int
]
lib.save_data_scalar_export0.restype = None

def save_data_scalar(filename: str, data: np.ndarray,
                     is_complex: bool, mesh, domain: np.ndarray,
                     w_points = None):
    """Save scalar field data to HDF5 file.

    Args:
        filename: Output filename
        data: nD array (complex or real)
        is_complex: Whether data is complex
        mesh: Mesh dimensions (e.g., [nx, ny, nz])
        domain: Domain vectors (2D array, shape [dimension, dimension])
        w_points: Frequency points (optional)
    """
    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    # Flatten and interleave data
    # For 3D arrays (k-space + frequency), use order='F' for w-k ordering
    # For 2D arrays, use default C order (row-major)
    order = 'F' if data.ndim == 3 else 'C'
    data_flat = data.ravel(order=order)
    if is_complex or np.iscomplexobj(data):
        data_interleaved = interleave_complex(data_flat)
        is_complex = True
    else:
        data_interleaved = np.asarray(data_flat, dtype=np.float32)

    total_size = data_flat.size
    mesh = np.asarray(mesh, dtype=np.int32)
    mesh_size = len(mesh)
    domain = np.asarray(domain, dtype=np.float32)
    domain_rows, domain_cols = domain.shape
    w_points = np.asarray(w_points, dtype=np.float32)
    w_size = len(w_points)

    # Flatten domain (row-major)
    domain_flat = domain.flatten()

    lib.save_data_scalar_export0(
        filename.encode('utf-8'),
        data_interleaved.ctypes.data_as(POINTER(c_float)),
        c_int(total_size), c_bool(is_complex),
        mesh.ctypes.data_as(POINTER(c_int)), c_int(mesh_size),
        domain_flat.ctypes.data_as(POINTER(c_float)),
        c_int(domain_rows), c_int(domain_cols),
        w_points.ctypes.data_as(POINTER(c_float)), c_int(w_size)
    )

lib.save_data_vector_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int
]
lib.save_data_vector_export0.restype = None

def save_data_vector(filename: str, data: np.ndarray,
                     nk: int, vec_len: int, is_complex: bool,
                     mesh, domain: np.ndarray,
                     w_points = None):
    """Save vector field data to HDF5 file.

    Args:
        filename: Output filename
        data: nD array (complex or real)
        nk: Number of k-points
        vec_len: Vector length (dimension)
        is_complex: Whether data is complex
        mesh: Mesh dimensions
        domain: Domain vectors
        w_points: Frequency points (optional)
    """
    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    # Flatten and interleave data
    # For 3D arrays (k-space + frequency), use order='F' for w-k ordering
    # For 2D arrays, use default C order (row-major)
    order = 'F' if data.ndim == 3 else 'C'
    data_flat = data.ravel(order=order)
    if is_complex or np.iscomplexobj(data):
        data_interleaved = interleave_complex(data_flat)
        is_complex = True
    else:
        data_interleaved = np.asarray(data_flat, dtype=np.float32)

    mesh = np.asarray(mesh, dtype=np.int32)
    mesh_size = len(mesh)
    domain = np.asarray(domain, dtype=np.float32)
    domain_rows, domain_cols = domain.shape
    w_points = np.asarray(w_points, dtype=np.float32)
    w_size = len(w_points)

    # Flatten domain
    domain_flat = domain.flatten()

    lib.save_data_vector_export0(
        filename.encode('utf-8'),
        data_interleaved.ctypes.data_as(POINTER(c_float)),
        c_int(nk), c_int(vec_len), c_bool(is_complex),
        mesh.ctypes.data_as(POINTER(c_int)), c_int(mesh_size),
        domain_flat.ctypes.data_as(POINTER(c_float)),
        c_int(domain_rows), c_int(domain_cols),
        w_points.ctypes.data_as(POINTER(c_float)), c_int(w_size)
    )

lib.save_data_matrix_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int
]
lib.save_data_matrix_export0.restype = None

def save_data_matrix(filename: str, data: np.ndarray,
                     num_states: int, mat_dim: int, is_complex: bool,
                     mesh, domain: np.ndarray,
                     w_points = None):
    """Save matrix field data to HDF5 file.

    Args:
        filename: Output filename
        data: nD array (complex or real)
        num_states: Number of states (matrix size)
        mat_dim: Matrix dimension
        is_complex: Whether data is complex
        mesh: Mesh dimensions
        domain: Domain vectors
        w_points: Frequency points (optional)
    """
    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    # Flatten and interleave data
    # For 3D arrays (k-space + frequency), use order='F' for w-k ordering
    # For 2D arrays, use default C order (row-major)
    order = 'F' if data.ndim == 3 else 'C'
    data_flat = data.ravel(order=order)
    if is_complex or np.iscomplexobj(data):
        data_interleaved = interleave_complex(data_flat)
        is_complex = True
    else:
        data_interleaved = np.asarray(data_flat, dtype=np.float32)

    mesh = np.asarray(mesh, dtype=np.int32)
    mesh_size = len(mesh)
    domain = np.asarray(domain, dtype=np.float32)
    domain_rows, domain_cols = domain.shape
    w_points = np.asarray(w_points, dtype=np.float32)
    w_size = len(w_points)

    # Flatten domain
    domain_flat = domain.flatten()
    num_matrices = data_flat.size // (mat_dim * mat_dim)

    lib.save_data_matrix_export0(
        filename.encode('utf-8'),
        data_interleaved.ctypes.data_as(POINTER(c_float)),
        c_int(num_matrices), c_int(mat_dim), c_bool(is_complex),
        mesh.ctypes.data_as(POINTER(c_int)), c_int(mesh_size),
        domain_flat.ctypes.data_as(POINTER(c_float)),
        c_int(domain_rows), c_int(domain_cols),
        w_points.ctypes.data_as(POINTER(c_float)), c_int(w_size)
    )

# BaseData exports
lib.BaseData_load.argtypes = [c_char_p]
lib.BaseData_load.restype = c_void_p

lib.BaseData_load_with_ordering.argtypes = [c_char_p, c_char_p]
lib.BaseData_load_with_ordering.restype = c_void_p

lib.BaseData_save.argtypes = [c_void_p, c_char_p]
lib.BaseData_save.restype = None

lib.BaseData_save_with_ordering.argtypes = [c_void_p, c_char_p, c_char_p]
lib.BaseData_save_with_ordering.restype = None

lib.destroy_BaseData.argtypes = [c_void_p]
lib.destroy_BaseData.restype = None

# BaseData metadata getters
lib.BaseData_get_is_complex.argtypes = [c_void_p]
lib.BaseData_get_is_complex.restype = c_int
lib.BaseData_get_is_vector.argtypes = [c_void_p]
lib.BaseData_get_is_vector.restype = c_int
lib.BaseData_get_is_matrix.argtypes = [c_void_p]
lib.BaseData_get_is_matrix.restype = c_int
lib.BaseData_get_with_k.argtypes = [c_void_p]
lib.BaseData_get_with_k.restype = c_int
lib.BaseData_get_with_w.argtypes = [c_void_p]
lib.BaseData_get_with_w.restype = c_int
lib.BaseData_get_as_mesh.argtypes = [c_void_p]
lib.BaseData_get_as_mesh.restype = c_int
lib.BaseData_get_n_indices.argtypes = [c_void_p]
lib.BaseData_get_n_indices.restype = c_int
lib.BaseData_get_dim_indices.argtypes = [c_void_p]
lib.BaseData_get_dim_indices.restype = c_int
lib.BaseData_get_dimension.argtypes = [c_void_p]
lib.BaseData_get_dimension.restype = c_int
lib.BaseData_get_nk.argtypes = [c_void_p]
lib.BaseData_get_nk.restype = c_int
lib.BaseData_get_nw.argtypes = [c_void_p]
lib.BaseData_get_nw.restype = c_int

# BaseData array getters
lib.BaseData_get_mesh_size.argtypes = [c_void_p]
lib.BaseData_get_mesh_size.restype = c_int
lib.BaseData_get_mesh.argtypes = [c_void_p, POINTER(c_int)]
lib.BaseData_get_mesh.restype = None
lib.BaseData_get_domain_rows.argtypes = [c_void_p]
lib.BaseData_get_domain_rows.restype = c_int
lib.BaseData_get_domain_cols.argtypes = [c_void_p]
lib.BaseData_get_domain_cols.restype = c_int
lib.BaseData_get_domain.argtypes = [c_void_p, POINTER(c_float)]
lib.BaseData_get_domain.restype = None
lib.BaseData_get_w_points_size.argtypes = [c_void_p]
lib.BaseData_get_w_points_size.restype = c_int
lib.BaseData_get_w_points.argtypes = [c_void_p, POINTER(c_float)]
lib.BaseData_get_w_points.restype = None

# BaseData data extraction
lib.BaseData_get_data_scalar.argtypes = [c_void_p, POINTER(c_float), POINTER(c_float)]
lib.BaseData_get_data_scalar.restype = None
lib.BaseData_get_data_matrix.argtypes = [c_void_p, POINTER(c_float), POINTER(c_float)]
lib.BaseData_get_data_matrix.restype = None

# Field get_data exports
lib.Field_R_get_data.argtypes = [c_void_p]
lib.Field_R_get_data.restype = c_void_p
lib.Field_C_get_data.argtypes = [c_void_p]
lib.Field_C_get_data.restype = c_void_p
lib.Field_RM_get_data.argtypes = [c_void_p]
lib.Field_RM_get_data.restype = c_void_p
lib.Field_CM_get_data.argtypes = [c_void_p]
lib.Field_CM_get_data.restype = c_void_p

class BaseData:
    """Python wrapper for BaseData C++ class with HDF5 save/load support."""

    def __init__(self, filename=None, ordering="k-w"):
        """
        Initialize BaseData from file or create empty instance.

        Args:
            filename: Path to HDF5 file to load (optional)
            ordering: Data ordering "k-w" or "w-k" (default: "k-w")
        """
        self.owns_ptr = True  # By default, we own the pointer
        if filename is None:
            self.ptr = None
            raise ValueError("BaseData requires a filename to load")
        else:
            if ordering == "k-w":
                self.ptr = lib.BaseData_load(c_char_p(filename.encode('utf-8')))
            else:
                self.ptr = lib.BaseData_load_with_ordering(
                    c_char_p(filename.encode('utf-8')),
                    c_char_p(ordering.encode('utf-8'))
                )

        if not self.ptr:
            raise RuntimeError('Failed to load BaseData')

        # Load metadata
        self._load_metadata()

    def _load_metadata(self):
        """Load metadata from C++ object."""
        self.is_complex = bool(lib.BaseData_get_is_complex(self.ptr))
        self.is_vector = bool(lib.BaseData_get_is_vector(self.ptr))
        self.is_matrix = bool(lib.BaseData_get_is_matrix(self.ptr))
        self.with_k = bool(lib.BaseData_get_with_k(self.ptr))
        self.with_w = bool(lib.BaseData_get_with_w(self.ptr))
        self.as_mesh = bool(lib.BaseData_get_as_mesh(self.ptr))
        self.n_indices = lib.BaseData_get_n_indices(self.ptr)
        self.dim_indices = lib.BaseData_get_dim_indices(self.ptr)
        self.dimension = lib.BaseData_get_dimension(self.ptr)
        self.nk = lib.BaseData_get_nk(self.ptr)
        self.nw = lib.BaseData_get_nw(self.ptr)

        # Load mesh
        mesh_size = lib.BaseData_get_mesh_size(self.ptr)
        if mesh_size > 0:
            mesh_buf = (c_int * mesh_size)()
            lib.BaseData_get_mesh(self.ptr, mesh_buf)
            self.mesh = np.array([mesh_buf[i] for i in range(mesh_size)], dtype=np.int32)
        else:
            self.mesh = np.array([], dtype=np.int32)

        # Load domain
        domain_rows = lib.BaseData_get_domain_rows(self.ptr)
        domain_cols = lib.BaseData_get_domain_cols(self.ptr)
        if domain_rows > 0 and domain_cols > 0:
            domain_buf = (c_float * (domain_rows * domain_cols))()
            lib.BaseData_get_domain(self.ptr, domain_buf)
            self.domain = np.array([domain_buf[i] for i in range(domain_rows * domain_cols)],
                                   dtype=np.float32).reshape(domain_rows, domain_cols)
        else:
            self.domain = np.array([], dtype=np.float32).reshape(0, 0)

        # Load w_points
        w_size = lib.BaseData_get_w_points_size(self.ptr)
        if w_size > 0:
            w_buf = (c_float * w_size)()
            lib.BaseData_get_w_points(self.ptr, w_buf)
            self.w_points = np.array([w_buf[i] for i in range(w_size)], dtype=np.float32)
        else:
            self.w_points = np.array([], dtype=np.float32)

    def save(self, filename, ordering="k-w"):
        """
        Save BaseData to HDF5 file.

        Args:
            filename: Output file path
            ordering: Data ordering "k-w" or "w-k" (default: "k-w")
        """
        if ordering == "k-w":
            lib.BaseData_save(self.ptr, c_char_p(filename.encode('utf-8')))
        else:
            lib.BaseData_save_with_ordering(
                self.ptr,
                c_char_p(filename.encode('utf-8')),
                c_char_p(ordering.encode('utf-8'))
            )

    def get_data(self):
        """
        Extract data as numpy array.

        Returns:
            numpy array with data (complex if is_complex=True)
        """
        if self.n_indices == 2:
            # Matrix data
            total_size = self.nk * self.nw * self.dim_indices * self.dim_indices
            real_buf = (c_float * total_size)()
            imag_buf = (c_float * total_size)() if self.is_complex else None

            lib.BaseData_get_data_matrix(self.ptr, real_buf, imag_buf if imag_buf else real_buf)

            real_data = np.array([real_buf[i] for i in range(total_size)], dtype=np.float32)
            if self.is_complex:
                imag_data = np.array([imag_buf[i] for i in range(total_size)], dtype=np.float32)
                data = real_data + 1j * imag_data
            else:
                data = real_data

            # Reshape to (nk*nw, dim_indices, dim_indices)
            return data.reshape(self.nk * self.nw, self.dim_indices, self.dim_indices)
        else:
            # Scalar data
            total_size = self.nk * self.nw
            real_buf = (c_float * total_size)()
            imag_buf = (c_float * total_size)() if self.is_complex else None

            lib.BaseData_get_data_scalar(self.ptr, real_buf, imag_buf if imag_buf else real_buf)

            real_data = np.array([real_buf[i] for i in range(total_size)], dtype=np.float32)
            if self.is_complex:
                imag_data = np.array([imag_buf[i] for i in range(total_size)], dtype=np.float32)
                data = real_data + 1j * imag_data
            else:
                data = real_data

            return data

    def __del__(self):
        try:
            if hasattr(self, 'ptr') and self.ptr and getattr(self, 'owns_ptr', True):
                lib.destroy_BaseData(self.ptr)
        except AttributeError:
            pass


# field = Field_R("sample_bands.dat")
# print(field(1, [-0.9, -0.9]))
#
# load_config("/home/g/Research/ffirefly/build/bin/input.cfg")
# print(epsilon(1, [0.1, 0.2, 0.3]))

class Field:
    """
    Field class that automatically dispatches to correct type based on field metadata.

    The field type is determined from the HDF5 file:
    - is_complex: True for complex fields, False for real
    - is_vector: True for vector fields (not currently used for matrix dispatch)
    - is_matrix: True for matrix fields, False for scalar

    Based on these flags, the appropriate return type is used when calling the field.
    """
    def __init__(self, filename=None):
        """
        Create Field from file or empty.
        
        Args:
            filename: Path to HDF5 field file (optional)
        """
        if filename is None:
            create = lib.Field_create
            create.restype = c_void_p
            self.ptr = create()
        else:
            create_from_file = lib.Field_from_file
            create_from_file.argtypes = [c_char_p]
            create_from_file.restype = c_void_p
            self.ptr = create_from_file(c_char_p(filename.encode('utf-8')))
        
        # Get type flags
        get_is_complex = lib.Field_is_complex
        get_is_complex.argtypes = [c_void_p]
        get_is_complex.restype = c_bool
        self.is_complex = get_is_complex(self.ptr)
        
        get_is_vector = lib.Field_is_vector
        get_is_vector.argtypes = [c_void_p]
        get_is_vector.restype = c_bool
        self.is_vector = get_is_vector(self.ptr)
        
        get_is_matrix = lib.Field_is_matrix
        get_is_matrix.argtypes = [c_void_p]
        get_is_matrix.restype = c_bool
        self.is_matrix = get_is_matrix(self.ptr)

        # Get plot metadata
        get_default_plot_type = lib.Field_get_default_plot_type
        get_default_plot_type.argtypes = [c_void_p]
        get_default_plot_type.restype = c_char_p
        self.default_plot_type = get_default_plot_type(self.ptr).decode('utf-8') if get_default_plot_type(self.ptr) else ""

        get_title = lib.Field_get_title
        get_title.argtypes = [c_void_p]
        get_title.restype = c_char_p
        self.title = get_title(self.ptr).decode('utf-8') if get_title(self.ptr) else ""

        get_x_label = lib.Field_get_x_label
        get_x_label.argtypes = [c_void_p]
        get_x_label.restype = c_char_p
        self.x_label = get_x_label(self.ptr).decode('utf-8') if get_x_label(self.ptr) else ""

        get_y_label = lib.Field_get_y_label
        get_y_label.argtypes = [c_void_p]
        get_y_label.restype = c_char_p
        self.y_label = get_y_label(self.ptr).decode('utf-8') if get_y_label(self.ptr) else ""
    
    def __del__(self):
        try:
            if hasattr(self, 'ptr') and self.ptr:
                destroy = lib.Field_destroy
                destroy.argtypes = [c_void_p]
                destroy(self.ptr)
                self.ptr = None
        except (AttributeError, OSError):
            pass
    
    def save(self, filename):
        """Save field to HDF5 file."""
        save_func = lib.Field_save
        save_func.argtypes = [c_void_p, c_char_p]
        save_func(self.ptr, c_char_p(filename.encode('utf-8')))
    
    def __call__(self, *args, **kwargs):
        """
        Evaluate field at point(s).
        
        Dispatches to appropriate method based on field type:
        - Scalar real: returns float or np.ndarray of floats
        - Scalar complex: returns complex or np.ndarray of complex
        - Matrix real: returns 2D np.ndarray of floats
        - Matrix complex: returns 2D np.ndarray of complex
        
        Args:
            k: k-point as array-like (x, y, z) or array of k-points
            w: frequency (optional, default 0.0)
        
        Or:
            w: frequency (for frequency-only evaluation)
        """
        if self.is_matrix and self.is_complex:
            return self._call_matrix_complex(*args, **kwargs)
        elif self.is_matrix and not self.is_complex:
            return self._call_matrix_real(*args, **kwargs)
        elif not self.is_matrix and self.is_complex:
            return self._call_scalar_complex(*args, **kwargs)
        else:
            return self._call_scalar_real(*args, **kwargs)
    
    def _call_scalar_real(self, k=None, w=0.0):
        """Call real scalar field."""
        if k is None:
            # Frequency-only evaluation
            call_func = lib.Field_call_scalar_real_w
            call_func.argtypes = [c_void_p, c_float]
            call_func.restype = c_float
            return call_func(self.ptr, c_float(w))
        
        k = np.atleast_2d(k).astype(np.float32)
        if k.shape[0] == 1:
            # Single point
            call_func = lib.Field_call_scalar_real_kw
            call_func.argtypes = [c_void_p, POINTER(c_float), c_int, c_float]
            call_func.restype = c_float
            k_ptr = k[0].ctypes.data_as(POINTER(c_float))
            return call_func(self.ptr, k_ptr, len(k[0]), c_float(w))
        else:
            # Multiple points
            call_func = lib.Field_call_scalar_real_list
            call_func.argtypes = [c_void_p, POINTER(c_float), c_int, c_int, c_float, POINTER(c_float)]
            k_flat = k.flatten()
            k_ptr = k_flat.ctypes.data_as(POINTER(c_float))
            result = np.zeros(k.shape[0], dtype=np.float32)
            result_ptr = result.ctypes.data_as(POINTER(c_float))
            call_func(self.ptr, k_ptr, k.shape[0], k.shape[1], c_float(w), result_ptr)
            return result
    
    def _call_scalar_complex(self, k=None, w=0.0):
        """Call complex scalar field."""
        if k is None:
            # Frequency-only evaluation
            call_func = lib.Field_call_scalar_complex_w
            call_func.argtypes = [c_void_p, c_float, POINTER(c_float), POINTER(c_float)]
            real_out = c_float()
            imag_out = c_float()
            call_func(self.ptr, c_float(w), byref(real_out), byref(imag_out))
            return complex(real_out.value, imag_out.value)
        
        k = np.atleast_2d(k).astype(np.float32)
        if k.shape[0] == 1:
            # Single point
            call_func = lib.Field_call_scalar_complex_kw
            call_func.argtypes = [c_void_p, POINTER(c_float), c_int, c_float, POINTER(c_float), POINTER(c_float)]
            k_ptr = k[0].ctypes.data_as(POINTER(c_float))
            real_out = c_float()
            imag_out = c_float()
            call_func(self.ptr, k_ptr, len(k[0]), c_float(w), byref(real_out), byref(imag_out))
            return complex(real_out.value, imag_out.value)
        else:
            # Multiple points
            call_func = lib.Field_call_scalar_complex_list
            call_func.argtypes = [c_void_p, POINTER(c_float), c_int, c_int, c_float, POINTER(c_float), POINTER(c_float)]
            k_flat = k.flatten()
            k_ptr = k_flat.ctypes.data_as(POINTER(c_float))
            real_out = np.zeros(k.shape[0], dtype=np.float32)
            imag_out = np.zeros(k.shape[0], dtype=np.float32)
            real_ptr = real_out.ctypes.data_as(POINTER(c_float))
            imag_ptr = imag_out.ctypes.data_as(POINTER(c_float))
            call_func(self.ptr, k_ptr, k.shape[0], k.shape[1], c_float(w), real_ptr, imag_ptr)
            return real_out + 1j * imag_out
    
    def _call_matrix_real(self, k, w=0.0):
        """Call real matrix field."""
        k = np.atleast_1d(k).astype(np.float32)
        call_func = lib.Field_call_matrix_real
        call_func.argtypes = [c_void_p, POINTER(c_float), c_int, c_float, POINTER(c_float), POINTER(c_int)]
        
        k_ptr = k.ctypes.data_as(POINTER(c_float))
        size_out = c_int()
        # Allocate maximum possible size
        max_size = 100
        result = np.zeros(max_size * max_size, dtype=np.float32)
        result_ptr = result.ctypes.data_as(POINTER(c_float))
        
        call_func(self.ptr, k_ptr, len(k), c_float(w), result_ptr, byref(size_out))
        
        n = size_out.value
        return result[:n*n].reshape(n, n)
    
    def _call_matrix_complex(self, k, w=0.0):
        """Call complex matrix field."""
        k = np.atleast_1d(k).astype(np.float32)
        call_func = lib.Field_call_matrix_complex
        call_func.argtypes = [c_void_p, POINTER(c_float), c_int, c_float, POINTER(c_float), POINTER(c_float), POINTER(c_int)]
        
        k_ptr = k.ctypes.data_as(POINTER(c_float))
        size_out = c_int()
        # Allocate maximum possible size
        max_size = 100
        real_out = np.zeros(max_size * max_size, dtype=np.float32)
        imag_out = np.zeros(max_size * max_size, dtype=np.float32)
        real_ptr = real_out.ctypes.data_as(POINTER(c_float))
        imag_ptr = imag_out.ctypes.data_as(POINTER(c_float))
        
        call_func(self.ptr, k_ptr, len(k), c_float(w), real_ptr, imag_ptr, byref(size_out))
        
        n = size_out.value
        return (real_out[:n*n] + 1j * imag_out[:n*n]).reshape(n, n)
    
    def get_data(self):
        """Get underlying BaseData object."""
        get_data_func = lib.Field_get_data
        get_data_func.argtypes = [c_void_p]
        get_data_func.restype = c_void_p
        ptr = get_data_func(self.ptr)
        return BaseData(ptr, owns_ptr=False)
