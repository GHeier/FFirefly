from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "fcode.fmodule",  # Name of the compiled submodule
        [
            "fcode/fmodule.cpp",
            "fcode/config/load/cpp_config.cpp",
            "fcode/config/load/c_config.c",
            "fcode/algorithms/interpolate.cpp",
            "fcode/objects/field.cpp",
            "fcode/objects/vec.cpp",
            "fcode/hamiltonian/potential.cpp",
            "fcode/hamiltonian/band_structure.cpp",
            "fcode/response/susceptibility.cpp",
            ], # Source file
        include_dirs=[
            'fcode', 
            'fcode/config', 
            'fcode/algorithms', 
            'fcode/objects', 
            'fcode/hamiltonian', 
            'fcode/response'
            ],  # Include directory for headers
        language="c++",  # Specify C++ compilation
    ),
]

setup(
    name="fcode",  # Name of the package
    version="0.1",
    packages=["fcode"],  # Include the source directory
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)
