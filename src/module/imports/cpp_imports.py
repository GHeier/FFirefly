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

    def get_faces_and_areas(self):
        """Returns tuple of (kpoints, areas) where kpoints is list of k-point vectors and areas is list of floats."""
        lib.Surface_num_faces_export0.argtypes = [ctypes.c_void_p]
        lib.Surface_num_faces_export0.restype = ctypes.c_int
        n = lib.Surface_num_faces_export0(self.ptr)

        if n <= 0:
            return ([], [])

        dims = (ctypes.c_int * n)()
        areas = (ctypes.c_float * n)()
        total_len = 3 * n  # Maximum possible size
        kpoints_buf = (ctypes.c_float * total_len)()

        lib.Surface_faces_and_areas_export0.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.Surface_faces_and_areas_export0.restype = None
        lib.Surface_faces_and_areas_export0(self.ptr, kpoints_buf, dims, areas, ctypes.c_int(n))

        # Reconstruct k-points list
        kpoints = []
        offset = 0
        for i in range(n):
            dim = dims[i]
            kpoints.append([kpoints_buf[offset + j] for j in range(dim)])
            offset += dim

        areas_list = [areas[i] for i in range(n)]
        return (kpoints, areas_list)

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
            if os.path.exists(filename):
                self.ptr = lib.Field_R_export2(c_char_p(filename.encode('utf-8')))
            else:
                raise FileNotFoundError(f"File {filename} does not exist")
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
        if len(bd.inds) == 2:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nk, bd.nw, bd.inds[0], bd.inds[1])
            return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        elif len(bd.inds) == 4:
            data_reshaped = data_array.reshape(bd.nk, bd.nw, bd.inds[0], bd.inds[1], bd.inds[2], bd.inds[3])
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
            if os.path.exists(filename):
                self.ptr = lib.Field_C_export2(c_char_p(filename.encode('utf-8')))
            else:
                raise FileNotFoundError(f"File {filename} does not exist")
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
        if len(bd.inds) == 2:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nk, bd.nw, bd.dim_indices, bd.dim_indices)
            return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        else:
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
lib.Hamiltonian_operator_export_list.argtypes = [c_void_p, POINTER(c_float), c_int, c_int, c_float, POINTER(c_float), POINTER(c_float), POINTER(c_int)]
lib.Hamiltonian_operator_export_list.restype = None
lib.Hamiltonian_file_found.argtypes = [c_void_p]
lib.Hamiltonian_file_found.restype = ctypes.c_bool

class Hamiltonian:
    def __init__(self):
        self.ptr = lib.Hamiltonian_export0()
        if not self.ptr:
            raise RuntimeError('Failed to initialize Hamiltonian')

    @property
    def file_found(self):
        """Check if Hamiltonian was loaded from file."""
        return lib.Hamiltonian_file_found(self.ptr)

    def __call__(self, k, w=0.0):
        """
        Evaluate Hamiltonian H(k, w) and return as numpy matrix or list of matrices.

        Args:
            k: momentum point (list or array) or list of momentum points
            w: frequency (default 0.0)

        Returns:
            numpy array of shape (n, n) with complex values, or list of such arrays
        """
        # Check if k is a list of points
        if isinstance(k, (list, tuple, np.ndarray)) and len(k) > 0:
            if isinstance(k[0], (list, tuple, np.ndarray)):
                # List of k-points
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

                # Allocate space for maximum size matrices
                max_size = 100
                real_output = (c_float * (num_points * max_size * max_size))()
                imag_output = (c_float * (num_points * max_size * max_size))()

                lib.Hamiltonian_operator_export_list(
                    self.ptr, points_flat, c_int(num_points), c_int(point_len), w_val,
                    real_output, imag_output, ctypes.byref(matrix_size)
                )

                n = matrix_size.value
                if n == 0:
                    return [np.array([[]], dtype=np.complex64) for _ in range(num_points)]

                # Reshape flattened arrays to list of matrices
                result = []
                for i in range(num_points):
                    offset = i * n * n
                    real_matrix = np.array([real_output[offset + j] for j in range(n*n)]).reshape(n, n)
                    imag_matrix = np.array([imag_output[offset + j] for j in range(n*n)]).reshape(n, n)
                    result.append(real_matrix + 1j * imag_matrix)

                return result

        # Single k-point
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

    def get_bands(self, k):
        """
        Diagonalize Hamiltonian at k-point(s) and return eigenvalues (band energies).

        Args:
            k: momentum point (list/array of length 3) or list of momentum points

        Returns:
            numpy array of eigenvalues (floats), or 2D array for multiple k-points
        """
        # Check if k is a list of points
        if isinstance(k, (list, tuple, np.ndarray)) and len(k) > 0:
            if isinstance(k[0], (list, tuple, np.ndarray)):
                # List of k-points
                points = k
                num_points = len(points)
                if num_points == 0:
                    return np.array([], dtype=np.float32)

                point_len = len(points[0])
                points_flat = (c_float * (num_points * point_len))()
                for i, p in enumerate(points):
                    for j, val in enumerate(p):
                        points_flat[i * point_len + j] = float(val)

                # Allocate space for eigenvalues
                max_bands = 100
                eigenvalues_out = (c_float * (num_points * max_bands))()
                num_bands = c_int(0)

                lib.Hamiltonian_get_bands_export_list(
                    self.ptr, points_flat, c_int(num_points), c_int(point_len),
                    eigenvalues_out, ctypes.byref(num_bands)
                )

                n = num_bands.value
                if n == 0:
                    return np.zeros((num_points, 0), dtype=np.float32)

                # Reshape to (num_points, n)
                result = np.array([eigenvalues_out[i] for i in range(num_points * n)], dtype=np.float32)
                return result.reshape(num_points, n)

        # Single k-point
        k_array = (c_float * len(k))(*[float(x) for x in k])
        k_len = c_int(len(k))

        # Allocate space for eigenvalues
        max_bands = 100
        eigenvalues_out = (c_float * max_bands)()
        num_bands = c_int(0)

        lib.Hamiltonian_get_bands_export0(
            self.ptr, k_array, k_len,
            eigenvalues_out, ctypes.byref(num_bands)
        )

        n = num_bands.value
        if n == 0:
            return np.array([], dtype=np.float32)

        return np.array([eigenvalues_out[i] for i in range(n)], dtype=np.float32)

    def get_wavefunctions(self, k):
        """
        Diagonalize Hamiltonian at k-point(s) and return eigenvalues and eigenvectors.

        Args:
            k: momentum point (list/array of length 3) or list of momentum points

        Returns:
            For single k-point:
                tuple (eigenvalues, eigenvectors) where:
                    - eigenvalues: numpy array of shape (n,) with floats
                    - eigenvectors: numpy array of shape (n, n) with complex values
                      (each column is an eigenvector)
            For multiple k-points:
                tuple (eigenvalues, eigenvectors) where:
                    - eigenvalues: numpy array of shape (num_points, n) with floats
                    - eigenvectors: numpy array of shape (num_points, n, n) with complex values
        """
        # Check if k is a list of points
        if isinstance(k, (list, tuple, np.ndarray)) and len(k) > 0:
            if isinstance(k[0], (list, tuple, np.ndarray)):
                # List of k-points
                points = k
                num_points = len(points)
                if num_points == 0:
                    return np.array([], dtype=np.float32), np.array([[]], dtype=np.complex64)

                point_len = len(points[0])
                points_flat = (c_float * (num_points * point_len))()
                for i, p in enumerate(points):
                    for j, val in enumerate(p):
                        points_flat[i * point_len + j] = float(val)

                # Allocate space for eigenvalues and eigenvectors
                max_bands = 100
                eigenvalues_out = (c_float * (num_points * max_bands))()
                eigvecs_real = (c_float * (num_points * max_bands * max_bands))()
                eigvecs_imag = (c_float * (num_points * max_bands * max_bands))()
                num_bands = c_int(0)

                lib.Hamiltonian_get_wavefunctions_export_list(
                    self.ptr, points_flat, c_int(num_points), c_int(point_len),
                    eigenvalues_out, eigvecs_real, eigvecs_imag,
                    ctypes.byref(num_bands)
                )

                n = num_bands.value
                if n == 0:
                    return np.zeros((num_points, 0), dtype=np.float32), np.zeros((num_points, 0, 0), dtype=np.complex64)

                # Extract eigenvalues
                eigs = np.array([eigenvalues_out[i] for i in range(num_points * n)], dtype=np.float32)
                eigs = eigs.reshape(num_points, n)

                # Extract eigenvectors (stored in column-major format for each k-point)
                vecs = np.zeros((num_points, n, n), dtype=np.complex64)
                for p in range(num_points):
                    offset = p * n * n
                    real_part = np.array([eigvecs_real[offset + i] for i in range(n*n)]).reshape(n, n)
                    imag_part = np.array([eigvecs_imag[offset + i] for i in range(n*n)]).reshape(n, n)
                    vecs[p] = real_part + 1j * imag_part

                return eigs, vecs

        # Single k-point
        k_array = (c_float * len(k))(*[float(x) for x in k])
        k_len = c_int(len(k))

        # Allocate space for eigenvalues and eigenvectors
        max_bands = 100
        eigenvalues_out = (c_float * max_bands)()
        eigvecs_real = (c_float * (max_bands * max_bands))()
        eigvecs_imag = (c_float * (max_bands * max_bands))()
        num_bands = c_int(0)

        lib.Hamiltonian_get_wavefunctions_export0(
            self.ptr, k_array, k_len,
            eigenvalues_out, eigvecs_real, eigvecs_imag,
            ctypes.byref(num_bands)
        )

        n = num_bands.value
        if n == 0:
            return np.array([], dtype=np.float32), np.array([[]], dtype=np.complex64)

        # Extract eigenvalues
        eigs = np.array([eigenvalues_out[i] for i in range(n)], dtype=np.float32)

        # Extract eigenvectors (stored in column-major format)
        real_part = np.array([eigvecs_real[i] for i in range(n*n)]).reshape(n, n)
        imag_part = np.array([eigvecs_imag[i] for i in range(n*n)]).reshape(n, n)
        vecs = real_part + 1j * imag_part

        return eigs, vecs

    def get_fermi_velocity(self, k):
        """
        Calculate Fermi velocity v_n(k) = ∇_k E_n(k) for all bands at k-point(s).

        Args:
            k: momentum point (list/array of length 2 or 3) or list of momentum points

        Returns:
            For single k-point:
                numpy array of shape (n, 3) where n is number of bands
                Each row is the velocity vector [vx, vy, vz] for that band
            For multiple k-points:
                numpy array of shape (num_points, n, 3)
        """
        # Check if k is a list of points
        if isinstance(k, (list, tuple, np.ndarray)) and len(k) > 0:
            if isinstance(k[0], (list, tuple, np.ndarray)):
                # List of k-points
                points = k
                num_points = len(points)
                if num_points == 0:
                    return np.zeros((0, 0, 3), dtype=np.float32)

                point_len = len(points[0])
                points_flat = (c_float * (num_points * point_len))()
                for i, p in enumerate(points):
                    for j, val in enumerate(p):
                        points_flat[i * point_len + j] = float(val)

                # Allocate space for velocities (num_points × max_bands × 3)
                max_bands = 100
                velocities_out = (c_float * (num_points * max_bands * 3))()
                num_bands = c_int(0)

                lib.Hamiltonian_get_fermi_velocity_export_list.argtypes = [
                    ctypes.c_void_p,
                    ctypes.POINTER(ctypes.c_float),
                    ctypes.c_int,
                    ctypes.c_int,
                    ctypes.POINTER(ctypes.c_float),
                    ctypes.POINTER(ctypes.c_int),
                ]
                lib.Hamiltonian_get_fermi_velocity_export_list.restype = None

                lib.Hamiltonian_get_fermi_velocity_export_list(
                    self.ptr, points_flat, c_int(num_points), c_int(point_len),
                    velocities_out, ctypes.byref(num_bands)
                )

                n = num_bands.value
                if n == 0:
                    return np.zeros((num_points, 0, 3), dtype=np.float32)

                # Reshape to (num_points, n, 3)
                result = np.zeros((num_points, n, 3), dtype=np.float32)
                for p in range(num_points):
                    for i in range(n):
                        idx = (p * n + i) * 3
                        result[p, i, 0] = velocities_out[idx + 0]
                        result[p, i, 1] = velocities_out[idx + 1]
                        result[p, i, 2] = velocities_out[idx + 2]

                return result

        # Single k-point
        k_array = (c_float * len(k))(*[float(x) for x in k])
        k_len = c_int(len(k))

        # Allocate space for velocities (max_bands × 3)
        max_bands = 100
        velocities_out = (c_float * (max_bands * 3))()
        num_bands = c_int(0)

        lib.Hamiltonian_get_fermi_velocity_export0.argtypes = [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.Hamiltonian_get_fermi_velocity_export0.restype = None

        lib.Hamiltonian_get_fermi_velocity_export0(
            self.ptr, k_array, k_len,
            velocities_out, ctypes.byref(num_bands)
        )

        n = num_bands.value
        if n == 0:
            return np.zeros((0, 3), dtype=np.float32)

        # Reshape to (n, 3)
        result = np.zeros((n, 3), dtype=np.float32)
        for i in range(n):
            result[i, 0] = velocities_out[i * 3 + 0]
            result[i, 1] = velocities_out[i * 3 + 1]
            result[i, 2] = velocities_out[i * 3 + 2]

        return result

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
lib.Field_RM_operator_export_w.argtypes = [c_void_p, c_float, POINTER(c_float), POINTER(c_int)]
lib.Field_RM_operator_export_w.restype = None
lib.Field_RM_operator_export_w_list.argtypes = [c_void_p, POINTER(c_float), c_int, POINTER(c_float), POINTER(c_int)]
lib.Field_RM_operator_export_w_list.restype = None

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
            if os.path.exists(filename):
                self.ptr = lib.Field_RM_export2(c_char_p(filename.encode('utf-8')))
            else:
                raise FileNotFoundError(f"File {filename} does not exist")
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

    def __call__(self, k=None, w=0.0):
        """
        Evaluate Field_RM at point(s) k with frequency w and return as numpy matrix or list of matrices.

        Args:
            k: momentum point (list or array), list of momentum points, or w value(s) if omitted
            w: frequency (default 0.0)

        Returns:
            numpy array of shape (n, n) with real values, or list of such arrays
        """
        # Overload for (w: float) - w-only call
        if k is not None and isinstance(k, (int, float, np.float32, np.float64)) and w == 0.0:
            w_val = c_float(k)
            matrix_size = c_int(0)

            # Allocate space for a maximum size matrix
            max_size = 100
            result = (c_float * (max_size * max_size))()

            lib.Field_RM_operator_export_w(
                self.ptr, w_val, result, ctypes.byref(matrix_size)
            )

            n = matrix_size.value
            if n == 0:
                return np.array([[]], dtype=np.float32)

            # Reshape flattened array to matrix
            matrix = np.array([result[i] for i in range(n*n)]).reshape(n, n)

            return matrix

        # Overload for (w_points: list[float]) - list of w values
        if k is not None and isinstance(k, (list, tuple, np.ndarray)) and len(k) > 0 and isinstance(k[0], (int, float, np.float32, np.float64)) and w == 0.0:
            w_points = k
            num_w = len(w_points)
            w_array = (c_float * num_w)(*[float(wval) for wval in w_points])
            matrix_size = c_int(0)

            # Allocate space for multiple matrices
            max_size = 100
            output = (c_float * (num_w * max_size * max_size))()

            lib.Field_RM_operator_export_w_list(
                self.ptr, w_array, c_int(num_w),
                output, ctypes.byref(matrix_size)
            )

            n = matrix_size.value
            if n == 0:
                return [np.array([[]], dtype=np.float32) for _ in range(num_w)]

            # Reshape to list of matrices
            matrices = []
            for i in range(num_w):
                start_idx = i * n * n
                end_idx = (i + 1) * n * n
                matrix = np.array([output[j] for j in range(start_idx, end_idx)]).reshape(n, n)
                matrices.append(matrix)

            return matrices

        # Check if k is a list of points (list of lists)
        if k is not None and len(k) > 0 and isinstance(k[0], (list, tuple, np.ndarray)):
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
lib.Field_CM_operator_export_w.argtypes = [c_void_p, c_float, POINTER(c_float), POINTER(c_float), POINTER(c_int)]
lib.Field_CM_operator_export_w.restype = None
lib.Field_CM_operator_export_w_list.argtypes = [c_void_p, POINTER(c_float), c_int, POINTER(c_float), POINTER(c_float), POINTER(c_int)]
lib.Field_CM_operator_export_w_list.restype = None

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
            if os.path.exists(filename):
                self.ptr = lib.Field_CM_export2(c_char_p(filename.encode('utf-8')))
            else:
                raise FileNotFoundError(f"File {filename} does not exist")
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

    def __call__(self, k=None, w=0.0):
        """
        Evaluate Field_CM at point(s) k with frequency w and return as numpy matrix or list of matrices.

        Args:
            k: momentum point (list or array), list of momentum points, or w value(s) if omitted
            w: frequency (default 0.0)

        Returns:
            numpy array of shape (n, n) with complex values, or list of such arrays
        """
        # Overload for (w: float) - w-only call
        if k is not None and isinstance(k, (int, float, np.float32, np.float64)) and w == 0.0:
            w_val = c_float(k)
            matrix_size = c_int(0)

            # Allocate space for a maximum size matrix
            max_size = 100
            real_result = (c_float * (max_size * max_size))()
            imag_result = (c_float * (max_size * max_size))()

            lib.Field_CM_operator_export_w(
                self.ptr, w_val, real_result, imag_result, ctypes.byref(matrix_size)
            )

            n = matrix_size.value
            if n == 0:
                return np.array([[]], dtype=np.complex64)

            # Reshape flattened arrays to matrix
            real_matrix = np.array([real_result[i] for i in range(n*n)]).reshape(n, n)
            imag_matrix = np.array([imag_result[i] for i in range(n*n)]).reshape(n, n)

            return real_matrix + 1j * imag_matrix

        # Overload for (w_points: list[float]) - list of w values
        if k is not None and isinstance(k, (list, tuple, np.ndarray)) and len(k) > 0 and isinstance(k[0], (int, float, np.float32, np.float64)) and w == 0.0:
            w_points = k
            num_w = len(w_points)
            w_array = (c_float * num_w)(*[float(wval) for wval in w_points])
            matrix_size = c_int(0)

            # Allocate space for multiple matrices
            max_size = 100
            real_output = (c_float * (num_w * max_size * max_size))()
            imag_output = (c_float * (num_w * max_size * max_size))()

            lib.Field_CM_operator_export_w_list(
                self.ptr, w_array, c_int(num_w),
                real_output, imag_output, ctypes.byref(matrix_size)
            )

            n = matrix_size.value
            if n == 0:
                return [np.array([[]], dtype=np.complex64) for _ in range(num_w)]

            # Reshape to list of matrices
            matrices = []
            for i in range(num_w):
                start_idx = i * n * n
                end_idx = (i + 1) * n * n
                real_mat = np.array([real_output[j] for j in range(start_idx, end_idx)]).reshape(n, n)
                imag_mat = np.array([imag_output[j] for j in range(start_idx, end_idx)]).reshape(n, n)
                matrices.append(real_mat + 1j * imag_mat)

            return matrices

        # Check if k is a list of points (list of lists)
        if k is not None and len(k) > 0 and isinstance(k[0], (list, tuple, np.ndarray)):
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
        if len(bd.inds) == 4:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nw, bd.nk, bd.inds[0], bd.inds[1], bd.inds[2], bd.inds[3])
            return data_reshaped
            #return np.moveaxis(data_reshaped, [0, 1], [1, 0])  # swap k and w axes
        if len(bd.inds) == 2:
            # Matrix: (nk*nw, dim, dim) -> (nw, nk, dim, dim)
            data_reshaped = data_array.reshape(bd.nw, bd.nk, bd.inds[0], bd.inds[1])
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
    if w_points is None:
        w_points = []

    with_w = len(w_points) > 0
    found_mesh = values.shape
    nbnd = 1

    if with_w and found_mesh[0] != len(w_points) or (not with_w and found_mesh[0] != mesh[0]):
        nbnd = found_mesh[0]

    domain = np.reshape(domain, -1).astype(np.float32)
    mesh_arr = np.array(mesh, dtype=np.int32)
    len_mesh = len(mesh_arr)

    is_complex = np.iscomplexobj(values)
    is_vector = False  # Vector support not implemented
    with_n = nbnd > 1

    if not is_complex:
        values_c = np.reshape(values, -1).astype(np.float32)
    else:
        values_c = interleave_complex(values)

    w_points = np.array(w_points, dtype=np.float32) if with_w else np.array([], dtype=np.float32)

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

# Set up the function signature
lib.load_config_export0.argtypes = [ctypes.c_char_p]
lib.load_config_export0.restype = None  # Equivalent to Cvoid

# Define the wrapper function
def load_config(path: str) -> None:
    lib.load_config_export0(path.encode("utf-8"))

# Save data functions
def save_data(filename: str, data: np.ndarray, mesh=None, domain=None,
              w_points=None, inds=None, n_indices=None, dim_indices=None):
    """Save data to HDF5 file with automatic dispatch based on data type.

    Args:
        filename: Output filename
        data: nD array (complex or real)
        mesh: Mesh dimensions (default: empty array)
        domain: Domain vectors (default: zeros)
        w_points: Frequency points (default: empty array)
        inds: Tensor index dimensions (e.g., [3,3] for 3x3 matrix)
        n_indices: DEPRECATED - Number of band indices (for backward compatibility)
        dim_indices: DEPRECATED - Dimension of indices (for backward compatibility)
    """
    # Set defaults
    if mesh is None:
        mesh = np.array([], dtype=np.int32)
    if domain is None:
        domain = np.zeros((0, 0), dtype=np.float32)
    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    is_complex = np.iscomplexobj(data)

    # Handle legacy API (n_indices, dim_indices) - convert to inds
    if inds is None and n_indices is not None:
        if n_indices == 0 or n_indices == 1:
            inds = []
        elif n_indices == 2:
            inds = [dim_indices if dim_indices else 1, dim_indices if dim_indices else 1]
        elif n_indices == 3:
            dim = dim_indices if dim_indices else 1
            inds = [dim, dim, dim]
        elif n_indices == 4:
            dim = dim_indices if dim_indices else 1
            inds = [dim, dim, dim, dim]
    elif inds is None and dim_indices is not None and dim_indices > 1:
        # Legacy: dim_indices > 1 implies matrix
        inds = [dim_indices, dim_indices]
    elif inds is None:
        # Default: scalar
        inds = []

    # Convert inds to list if needed
    if isinstance(inds, np.ndarray):
        inds = inds.tolist()

    rank = len(inds)

    # Determine data type based on rank
    if rank == 4:
        # 4D tensor data
        save_data_tensor4(filename, data, is_complex, mesh, domain, w_points, inds)
    elif rank == 3:
        # 3D tensor data
        save_data_tensor3(filename, data, is_complex, mesh, domain, w_points, inds)
    elif rank == 2:
        # Matrix data
        save_data_matrix(filename, data, is_complex, mesh, domain, w_points, inds)
    elif rank == 1:
        # Vector data
        save_data_vector(filename, data, is_complex, mesh, domain, w_points, inds)
    else:
        # Scalar data (rank == 0)
        save_data_scalar(filename, data, is_complex, mesh, domain, w_points)

lib.save_data_scalar_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int,
    POINTER(c_float), c_int, c_int
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
    # Use C order (row-major) for w-k ordering where w varies slowest
    # For (nw, nky, nkx) array, C order gives correct w-k ordering
    data_flat = data.ravel(order='C')
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
        w_points.ctypes.data_as(POINTER(c_float)), c_int(w_size),
        None, c_int(0), c_int(0)
    )

lib.save_data_vector_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int,
    POINTER(c_float), c_int, c_int
]
lib.save_data_vector_export0.restype = None

def save_data_vector(filename: str, data: np.ndarray,
                     nk_or_is_complex = None, vec_len_or_mesh = None,
                     is_complex_or_domain = None, mesh_or_w_points = None,
                     domain_or_inds = None, w_points = None, inds = None):
    """Save vector field data to HDF5 file.

    Supports two API styles:
    - New API: save_data_vector(filename, data, is_complex, mesh, domain, w_points=None, inds=None)
    - Old API: save_data_vector(filename, data, nk, vec_len, is_complex, mesh, domain, w_points=None)

    Args:
        filename: Output filename
        data: nD array (complex or real)
        For new API: is_complex, mesh, domain, w_points, inds
        For old API: nk, vec_len, is_complex, mesh, domain, w_points
    """
    # Detect which API is being used based on parameter types
    if isinstance(nk_or_is_complex, bool):
        # New API: (filename, data, is_complex, mesh, domain, w_points, inds)
        is_complex = nk_or_is_complex
        mesh = vec_len_or_mesh
        domain = is_complex_or_domain
        w_points = mesh_or_w_points
        inds = domain_or_inds
    elif isinstance(nk_or_is_complex, int):
        # Old API: (filename, data, nk, vec_len, is_complex, mesh, domain, w_points)
        nk = nk_or_is_complex
        vec_len = vec_len_or_mesh
        is_complex = is_complex_or_domain
        mesh = mesh_or_w_points
        domain = domain_or_inds
        # w_points is the 8th parameter (already set)
        # inds will be derived from vec_len
        inds = [vec_len]
    else:
        raise ValueError("Invalid arguments to save_data_vector")

    # Use same save mechanism as scalar - vectors are saved as scalars
    save_data_scalar(filename, data, is_complex, mesh, domain, w_points)

lib.save_data_matrix_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int,
    POINTER(c_float), c_int, c_int
]
lib.save_data_matrix_export0.restype = None

lib.save_data_tensor3_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int,
    POINTER(c_float), c_int, c_int
]
lib.save_data_tensor3_export0.restype = None

lib.save_data_tensor4_export0.argtypes = [
    c_char_p, POINTER(c_float), c_int, c_int, c_bool,
    POINTER(c_int), c_int, POINTER(c_float), c_int, c_int, POINTER(c_float), c_int,
    POINTER(c_float), c_int, c_int
]
lib.save_data_tensor4_export0.restype = None

def save_data_matrix(filename: str, data: np.ndarray,
                     num_matrices_or_is_complex = None, mat_dim_or_mesh = None,
                     is_complex_or_domain = None, mesh_or_w_points = None,
                     domain_or_inds = None, w_points = None, inds = None):
    """Save matrix field data to HDF5 file.

    Supports two API styles:
    - New API: save_data_matrix(filename, data, is_complex, mesh, domain, w_points=None, inds=None)
    - Old API: save_data_matrix(filename, data, num_matrices, mat_dim, is_complex, mesh, domain, w_points=None)

    Args:
        filename: Output filename
        data: nD array (complex or real)
        For new API: is_complex, mesh, domain, w_points, inds
        For old API: num_matrices, mat_dim, is_complex, mesh, domain, w_points
    """
    # Detect which API is being used based on parameter types
    if isinstance(num_matrices_or_is_complex, bool):
        # New API: (filename, data, is_complex, mesh, domain, w_points, inds)
        is_complex = num_matrices_or_is_complex
        mesh = mat_dim_or_mesh
        domain = is_complex_or_domain
        w_points = mesh_or_w_points
        inds = domain_or_inds
        mat_dim = None  # Will be derived from inds
        num_matrices = None  # Will be derived from data size
    elif isinstance(num_matrices_or_is_complex, int):
        # Old API: (filename, data, num_matrices, mat_dim, is_complex, mesh, domain, w_points)
        num_matrices = num_matrices_or_is_complex
        mat_dim = mat_dim_or_mesh
        is_complex = is_complex_or_domain
        mesh = mesh_or_w_points
        domain = domain_or_inds
        # w_points is the 8th parameter (already set)
        # inds will be derived from mat_dim
        inds = [mat_dim, mat_dim]
    else:
        raise ValueError("Invalid arguments to save_data_matrix")

    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    # Flatten and interleave data
    # Use C order (row-major) for w-k ordering where w varies slowest
    data_flat = data.ravel(order='C')
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

    # Derive mat_dim and num_matrices if not provided
    if mat_dim is None:
        mat_dim = inds[0]  # Assuming square matrices (inds[0] == inds[1])
    if num_matrices is None:
        num_matrices = data_flat.size // (mat_dim * mat_dim)

    lib.save_data_matrix_export0(
        filename.encode('utf-8'),
        data_interleaved.ctypes.data_as(POINTER(c_float)),
        c_int(num_matrices), c_int(mat_dim), c_bool(is_complex),
        mesh.ctypes.data_as(POINTER(c_int)), c_int(mesh_size),
        domain_flat.ctypes.data_as(POINTER(c_float)),
        c_int(domain_rows), c_int(domain_cols),
        w_points.ctypes.data_as(POINTER(c_float)), c_int(w_size),
        None, c_int(0), c_int(0)
    )

def save_data_tensor3(filename: str, data: np.ndarray,
                     is_complex: bool, mesh, domain: np.ndarray,
                     w_points = None, inds = None):
    """Save 3D tensor field data to HDF5 file.

    Args:
        filename: Output filename
        data: nD array (complex or real)
        is_complex: Whether data is complex
        mesh: Mesh dimensions
        domain: Domain vectors
        w_points: Frequency points (optional)
        inds: Tensor index dimensions (e.g., [3,3,3] for 3x3x3 tensor)
    """
    # Extract dimensions from inds
    if inds is None or len(inds) != 3:
        raise ValueError("save_data_tensor3 requires inds with 3 dimensions")

    # For now, assume all dimensions are equal (as C++ export expects)
    ten_dim = inds[0]
    if not all(d == ten_dim for d in inds):
        raise ValueError("save_data_tensor3 currently requires all tensor dimensions to be equal")

    # Calculate number of tensors from data shape
    tensor_size = ten_dim * ten_dim * ten_dim
    num_tensors = data.size // tensor_size
    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    # Flatten and interleave data
    # Use C order (row-major) for w-k ordering where w varies slowest
    data_flat = data.ravel(order='C')
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

    lib.save_data_tensor3_export0(
        filename.encode('utf-8'),
        data_interleaved.ctypes.data_as(POINTER(c_float)),
        c_int(num_tensors), c_int(ten_dim), c_bool(is_complex),
        mesh.ctypes.data_as(POINTER(c_int)), c_int(mesh_size),
        domain_flat.ctypes.data_as(POINTER(c_float)),
        c_int(domain_rows), c_int(domain_cols),
        w_points.ctypes.data_as(POINTER(c_float)), c_int(w_size),
        None, c_int(0), c_int(0)
    )

def save_data_tensor4(filename: str, data: np.ndarray,
                     is_complex: bool, mesh, domain: np.ndarray,
                     w_points = None, inds = None):
    """Save 4D tensor field data to HDF5 file.

    Args:
        filename: Output filename
        data: nD array (complex or real)
        is_complex: Whether data is complex
        mesh: Mesh dimensions
        domain: Domain vectors
        w_points: Frequency points (optional)
        inds: Tensor index dimensions (e.g., [3,3,3,3] for 3x3x3x3 tensor)
    """
    # Extract dimensions from inds
    if inds is None or len(inds) != 4:
        raise ValueError("save_data_tensor4 requires inds with 4 dimensions")

    # For now, assume all dimensions are equal (as C++ export expects)
    ten_dim = inds[0]
    if not all(d == ten_dim for d in inds):
        raise ValueError("save_data_tensor4 currently requires all tensor dimensions to be equal")

    # Calculate number of tensors from data shape
    tensor_size = ten_dim * ten_dim * ten_dim * ten_dim
    num_tensors = data.size // tensor_size
    if w_points is None:
        w_points = np.array([], dtype=np.float32)

    # Flatten and interleave data
    # Use C order (row-major) for w-k ordering where w varies slowest
    data_flat = data.ravel(order='C')
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

    lib.save_data_tensor4_export0(
        filename.encode('utf-8'),
        data_interleaved.ctypes.data_as(POINTER(c_float)),
        c_int(num_tensors), c_int(ten_dim), c_bool(is_complex),
        mesh.ctypes.data_as(POINTER(c_int)), c_int(mesh_size),
        domain_flat.ctypes.data_as(POINTER(c_float)),
        c_int(domain_rows), c_int(domain_cols),
        w_points.ctypes.data_as(POINTER(c_float)), c_int(w_size),
        None, c_int(0), c_int(0)
    )

# BaseData exports
lib.BaseData_load.argtypes = [c_char_p]
lib.BaseData_load.restype = c_void_p

lib.BaseData_load_with_ordering.argtypes = [c_char_p, c_char_p]
lib.BaseData_load_with_ordering.restype = c_void_p

lib.BaseData_save.argtypes = [c_void_p, c_char_p]
lib.BaseData_save.restype = None

# lib.BaseData_save_with_ordering.argtypes = [c_void_p, c_char_p, c_char_p]
# lib.BaseData_save_with_ordering.restype = None

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
lib.BaseData_get_rank.argtypes = [c_void_p]
lib.BaseData_get_rank.restype = c_int
lib.BaseData_get_inds_size.argtypes = [c_void_p]
lib.BaseData_get_inds_size.restype = c_int
lib.BaseData_get_inds.argtypes = [c_void_p, POINTER(c_int)]
lib.BaseData_get_inds.restype = None
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
lib.BaseData_get_points_rows.argtypes = [c_void_p]
lib.BaseData_get_points_rows.restype = c_int
lib.BaseData_get_points_cols.argtypes = [c_void_p]
lib.BaseData_get_points_cols.restype = c_int
lib.BaseData_get_points.argtypes = [c_void_p, POINTER(c_float)]
lib.BaseData_get_points.restype = None

# BaseData data extraction
lib.BaseData_get_data_scalar.argtypes = [c_void_p, POINTER(c_float), POINTER(c_float)]
lib.BaseData_get_data_scalar.restype = None
lib.BaseData_get_data_matrix.argtypes = [c_void_p, POINTER(c_float), POINTER(c_float)]
lib.BaseData_get_data_matrix.restype = None
lib.BaseData_get_data_tensor3.argtypes = [c_void_p, POINTER(c_float), POINTER(c_float)]
lib.BaseData_get_data_tensor3.restype = None
lib.BaseData_get_data_tensor4.argtypes = [c_void_p, POINTER(c_float), POINTER(c_float)]
lib.BaseData_get_data_tensor4.restype = None

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

    def __init__(self, filename=None, ordering="k-w", owns_ptr=True):
        """
        Initialize BaseData from file or from existing pointer.

        Args:
            filename: Path to HDF5 file to load, or an existing c_void_p pointer
            ordering: Data ordering "k-w" or "w-k" (default: "k-w")
            owns_ptr: Whether this object owns the pointer (for memory management)
        """
        self.owns_ptr = owns_ptr

        # Check if filename is actually a pointer (from Field.get_data())
        if isinstance(filename, int) or (hasattr(filename, 'value') and isinstance(filename, c_void_p)):
            # filename is actually a pointer
            self.ptr = filename
        elif filename is None:
            self.ptr = None
            raise ValueError("BaseData requires a filename to load")
        else:
            # filename is a string path
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

        # Load inds array
        inds_size = lib.BaseData_get_inds_size(self.ptr)
        if inds_size > 0:
            inds_buf = (c_int * inds_size)()
            lib.BaseData_get_inds(self.ptr, inds_buf)
            self.inds = list(inds_buf)
        else:
            self.inds = []

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

        # Load points
        points_rows = lib.BaseData_get_points_rows(self.ptr)
        points_cols = lib.BaseData_get_points_cols(self.ptr)
        if points_rows > 0 and points_cols > 0:
            points_buf = (c_float * (points_rows * points_cols))()
            lib.BaseData_get_points(self.ptr, points_buf)
            self.points = np.array([points_buf[i] for i in range(points_rows * points_cols)],
                                   dtype=np.float32).reshape(points_rows, points_cols)
        else:
            self.points = np.array([], dtype=np.float32).reshape(0, 0)

    def save(self, filename, ordering="k-w"):
        """
        Save BaseData to HDF5 file.

        Args:
            filename: Output file path
            ordering: Data ordering "k-w" or "w-k" (default: "k-w")
        """
        # Get data (flattened) and reshape for save_data
        data = self.data

        # Reshape data for save_data (expects w-k ordering with explicit dimensions)
        # For scalar: shape should be (nw, nk_total) where nk_total = product of mesh
        rank = len(self.inds)

        if rank == 0:
            # Scalar: reshape to (nw, *mesh_dims)
            mesh_shape = tuple(self.mesh) if len(self.mesh) > 0 else (self.nk,)
            shaped_data = data.reshape(self.nw, *mesh_shape)
        else:
            # Tensor: check if data is already properly shaped
            mesh_shape = tuple(self.mesh) if len(self.mesh) > 0 else (self.nk,)
            tensor_shape = tuple(self.inds)
            expected_shape = (self.nw, *mesh_shape, *tensor_shape)

            if data.shape == expected_shape:
                # Data is already properly shaped, use as-is
                shaped_data = data
            else:
                # Data needs reshaping (e.g., from flat array)
                shaped_data = data.reshape(self.nw, *mesh_shape, *tensor_shape)

        # Use save_data() to save the data
        save_data(filename, shaped_data, mesh=self.mesh, domain=self.domain,
                  w_points=self.w_points, inds=list(self.inds))

    @property
    def data(self):
        """
        Access data as numpy array.

        Returns:
            numpy array with data (complex if is_complex=True)
        """
        return self.get_data()

    @data.setter
    def data(self, new_data):
        """
        Set data array (updates internal C++ data).

        Args:
            new_data: numpy array with new data
        """
        self.set_data(new_data)

    def get_data(self):
        """
        Extract data as numpy array.

        Returns:
            numpy array with data (complex if is_complex=True)
        """
        # If data was set via setter, return that
        if hasattr(self, '_data'):
            return self._data

        rank = len(self.inds)

        # Calculate total tensor size
        tensor_size = 1
        for d in self.inds:
            tensor_size *= d

        total_size = self.nk * self.nw * tensor_size

        if rank == 4:
            # 4D tensor data
            real_buf = (c_float * total_size)()
            imag_buf = (c_float * total_size)() if self.is_complex else None

            lib.BaseData_get_data_tensor4(self.ptr, real_buf, imag_buf if imag_buf else real_buf)

            real_data = np.array([real_buf[i] for i in range(total_size)], dtype=np.float32)
            if self.is_complex:
                imag_data = np.array([imag_buf[i] for i in range(total_size)], dtype=np.float32)
                data = real_data + 1j * imag_data
            else:
                data = real_data

            # Reshape to (nk*nw, inds[0], inds[1], inds[2], inds[3])
            return data.reshape(self.nk * self.nw, self.inds[0], self.inds[1], self.inds[2], self.inds[3])
        elif rank == 3:
            # 3D tensor data
            real_buf = (c_float * total_size)()
            imag_buf = (c_float * total_size)() if self.is_complex else None

            lib.BaseData_get_data_tensor3(self.ptr, real_buf, imag_buf if imag_buf else real_buf)

            real_data = np.array([real_buf[i] for i in range(total_size)], dtype=np.float32)
            if self.is_complex:
                imag_data = np.array([imag_buf[i] for i in range(total_size)], dtype=np.float32)
                data = real_data + 1j * imag_data
            else:
                data = real_data

            # Reshape to (nk*nw, inds[0], inds[1], inds[2])
            return data.reshape(self.nk * self.nw, self.inds[0], self.inds[1], self.inds[2])
        elif rank == 2:
            # Matrix data
            real_buf = (c_float * total_size)()
            imag_buf = (c_float * total_size)() if self.is_complex else None

            lib.BaseData_get_data_matrix(self.ptr, real_buf, imag_buf if imag_buf else real_buf)

            real_data = np.array([real_buf[i] for i in range(total_size)], dtype=np.float32)
            if self.is_complex:
                imag_data = np.array([imag_buf[i] for i in range(total_size)], dtype=np.float32)
                data = real_data + 1j * imag_data
            else:
                data = real_data

            # Reshape to (nk*nw, inds[0], inds[1])
            return data.reshape(self.nk * self.nw, self.inds[0], self.inds[1])
        elif rank == 1:
            # Vector data - need to add export function for this
            real_buf = (c_float * total_size)()
            imag_buf = (c_float * total_size)() if self.is_complex else None

            # For now, use scalar export and reshape
            lib.BaseData_get_data_scalar(self.ptr, real_buf, imag_buf if imag_buf else real_buf)

            real_data = np.array([real_buf[i] for i in range(total_size)], dtype=np.float32)
            if self.is_complex:
                imag_data = np.array([imag_buf[i] for i in range(total_size)], dtype=np.float32)
                data = real_data + 1j * imag_data
            else:
                data = real_data

            # Reshape to (nk*nw, inds[0])
            return data.reshape(self.nk * self.nw, self.inds[0])
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

    def set_data(self, new_data):
        """
        Set data array (stored as Python attribute for use in save()).

        Args:
            new_data: numpy array with new data
        """
        new_data = np.asarray(new_data)

        # Validate shape matches BaseData properties
        rank = len(self.inds)
        tensor_size = 1
        for d in self.inds:
            tensor_size *= d

        expected_total = self.nk * self.nw * tensor_size
        if new_data.size != expected_total:
            raise ValueError(f"Data size mismatch: expected {expected_total}, got {new_data.size}")

        # Update is_complex based on new data type
        self.is_complex = np.iscomplexobj(new_data)

        # Flatten and store data as internal attribute (for consistency with get_data)
        self._data = new_data.ravel()

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


def get_reduced_grid(grid, lattice="SC"):
    """
    Get reduced k-point grid using symmetry equivalence classes.

    Parameters:
    -----------
    grid : list of int
        Grid dimensions, e.g., [5, 5] for 2D or [5, 5, 5] for 3D
    lattice : str
        Lattice type, e.g., "SC" (simple cubic), "BCC", "FCC"

    Returns:
    --------
    list of list of list of int
        Nested list structure: reduced_grid[group][point][coordinate]
        Each group contains symmetry-equivalent k-points
    """
    grid_size = len(grid)
    grid_arr = (c_int * grid_size)(*grid)
    lattice_cstr = lattice.encode('utf-8')

    # Allocate output buffers (maximum possible size)
    prod = 1
    for g in grid:
        prod *= g

    indices_out = (c_int * (prod * 3))()  # max 3 dimensions per point
    group_sizes = (c_int * prod)()  # max prod groups (worst case: no symmetry)
    point_dims = (c_int * prod)()  # dimension for each point
    num_groups = c_int()
    total_points = c_int()

    # Set up function
    lib.get_reduced_grid_export0.argtypes = [
        POINTER(c_int), c_int,  # grid, grid_size
        c_char_p,  # lattice
        POINTER(c_int),  # indices_out
        POINTER(c_int),  # group_sizes
        POINTER(c_int),  # point_dims
        POINTER(c_int),  # num_groups
        POINTER(c_int)   # total_points
    ]
    lib.get_reduced_grid_export0.restype = None

    # Call C++ function
    lib.get_reduced_grid_export0(
        grid_arr, grid_size,
        lattice_cstr,
        indices_out, group_sizes, point_dims,
        byref(num_groups), byref(total_points)
    )

    # Reconstruct nested structure
    result = []
    offset = 0
    group_start = 0

    for i in range(num_groups.value):
        group = []
        n_points = group_sizes[i]

        for j in range(n_points):
            point_idx = group_start + j
            dim = point_dims[point_idx]
            point = [indices_out[point_idx * 3 + k] for k in range(dim)]
            group.append(point)

        result.append(group)
        group_start += n_points

    return result
