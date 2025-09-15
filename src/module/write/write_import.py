from . import generate_python_wrappers
from . import generate_julia_wrappers
from . import base
import copy


def write_py_import(funcs, classes):
    with open("imports/cpp_imports.py", "r") as file:
        lines = file.readlines()
    newlines = generate_python_wrappers.create_python_function(funcs, classes)
    lines = base.remove_lines_between_phrases(
        lines, "# Begin Functions", "# End Functions"
    )
    lines = base.add_lines_between_phrases(
        lines, newlines, "# Begin Functions", "# End Functions"
    )
    with open("imports/cpp_imports.py", "w") as file:
        file.writelines(lines)


def write_jl_import(INPUT):
    input = copy.copy(INPUT)
    with open("imports/cpp_imports.jl", "r") as file:
        lines = file.readlines()
    newlines = generate_julia_wrappers.create_julia_function(INPUT)
    lines = base.remove_lines_between_phrases(
        lines, "# Begin Functions", "# End Functions"
    )
    lines = base.add_lines_between_phrases(
        lines, newlines, "# Begin Functions", "# End Functions"
    )
    with open("imports/cpp_imports.jl", "w") as file:
        file.writelines(lines)


def write_import(funcs, classes):
    write_py_import(funcs, classes)
    # write_jl_import(INPUT)
