from . import generate_wrappers
from . import base
import copy


def write_import(INPUT):
    input = copy.copy(INPUT)
    with open("imports/cpp_imports.py", "r") as file:
        lines = file.readlines()
    newlines = generate_wrappers.create_python_function(INPUT)
    lines = base.remove_lines_between_phrases(
        lines, "# Begin Functions", "# End Functions"
    )
    lines = base.add_lines_between_phrases(
        lines, newlines, "# Begin Functions", "# End Functions"
    )
    with open("imports/cpp_imports.py", "w") as file:
        file.writelines(lines)
    # generate_wrappers.create_function(INPUT)
