from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import pybind11
import os

# Helper to find source files for the pybind11 module
def collect_cpp_sources(base_dir, extensions=[".cpp", ".c"]):
    sources = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(tuple(extensions)):
                sources.append(os.path.join(root, file))
    return sources

# Collect sources for the C++ module
cpp_sources = ["src/config/load/pybind_module.cpp"]

# Define the pybind11 extension module
ext_modules = [
    Extension(
        "fcode",  # Name of the pybind11 module
        cpp_sources,  # C++ source files
        include_dirs=[
            pybind11.get_include(),
            "src/hamiltonian",  # Path to Hamiltonian headers
            "src/algorithms",  # Path to Algorithms headers
        ],
        language="c++",
    )
]

# Setup script
setup(
    name="fcode",
    version="0.1",
    description="A hybrid Python package with C++ bindings and Python scripts",
    packages=find_packages(where="src"),  # Automatically detect packages in src/
    package_dir={"": "src"},  # Specify src/ as the root for packages
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
