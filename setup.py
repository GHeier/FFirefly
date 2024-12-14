from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "fcode.fmodule",  # Name of the compiled submodule
        ["fcode/fmodule.cpp"], # Source file
        include_dirs=['fcode'],  # Include directory for headers
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
