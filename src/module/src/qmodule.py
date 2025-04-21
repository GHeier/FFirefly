from ctypes import CDLL, c_void_p, c_char_p, c_float, POINTER, Structure, byref
import numpy as np

class Field_C_Wrapper:
    def __init__(self):
        # Load the shared library
        self.lib = CDLL("/home/g/Research/ffirefly/build/lib/libfcode.so")
        
        # Define function prototypes
        self._setup_function_prototypes()
        
    def _setup_function_prototypes(self):
        # Field_C creation/destruction
        self.lib.create_Field_C.restype = c_void_p
        self.lib.load_Field_C.argtypes = [c_char_p]
        self.lib.load_Field_C.restype = c_void_p
        self.lib.destroy_Field_C.argtypes = [c_void_p]

        # Field_C operator() calls
        self.lib.Field_C_call.argtypes = [
            c_void_p,           # Field_C pointer
            POINTER(c_float),    # point array
            c_float,             # w
            c_int,               # len
            POINTER(c_float),    # real_result
            POINTER(c_float)     # imag_result
        ]

        self.lib.Field_C_call_w.argtypes = [
            c_void_p,
            c_float,
            POINTER(c_float),
            POINTER(c_float)
        ]

    class Field_C:
        def __init__(self, ptr):
            self.ptr = c_void_p(ptr)
            
        def __call__(self, point=None, w=0.0):
            if point is not None:
                return self._call_with_point(point, w)
            else:
                return self._call_with_w(w)
                
        def _call_with_point(self, point, w):
            real_result = c_float()
            imag_result = c_float()
            arr = (c_float * len(point))(*point)
            
            self.lib.Field_C_call(
                self.ptr,
                arr,
                c_float(w),
                len(point),
                byref(real_result),
                byref(imag_result)
            )
            return complex(real_result.value, imag_result.value)
            
        def _call_with_w(self, w):
            real_result = c_float()
            imag_result = c_float()
            self.lib.Field_C_call_w(
                self.ptr,
                c_float(w),
                byref(real_result),
                byref(imag_result))
            return complex(real_result.value, imag_result.value)
            
        def __del__(self):
            self.lib.destroy_Field_C(self.ptr)

    # Factory methods
    def create_field(self):
        return self.Field_C(self.lib.create_Field_C())
    
    def load_field(self, filename):
        return self.Field_C(self.lib.load_Field_C(filename.encode()))

field = Field_C_Wrapper()

# Create a Field_C object
f = field.load_field("/home/g/Research/Materials/Test/test_chi.dat")
