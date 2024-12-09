from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "hamiltonian.fcode",  # Fully qualified module name
        ["hamiltonian/band_structure.cpp"],  # Path to the C++ source file
    ),
]

setup(
    name="fcode",
    version="0.1",
    packages=["src"],  # Include all Python modules under "src"
    ext_modules=ext_modules,  # Include the C++ extension
    cmdclass={"build_ext": build_ext},
)

