from . import base
import copy
from collections import defaultdict

letters = "abcdefghijklmnopq"
variables = {"Vec": "v"}


def arg_cpp_to_c(arg, name):
    if arg == "vector" or arg == "Vec":
        return ["float*", "int"]
    if (
        arg == "vector"
        or arg == "Vec"
        or arg == "vector<Vec>"
        or arg == "vector<vector<float>>"
    ):
        return ""
    if arg == "string":
        return "char*"
    if arg == "ptr":
        return name + "*"
    return arg


def type_convert_lines(a, i):
    lines = []
    if a == "Vec":
        lines.append(
            "    Vec "
            + variables["Vec"]
            + "("
            + letters[i - 1]
            + ", "
            + letters[i]
            + ");"
        )
    return lines


def filter_default_args(args):
    arguments = []
    for x in args:
        arguments.append(x.split("=")[0])
    return arguments


def create_declaration(name, args, ret, num):
    arguments = filter_default_args(args)

    line = 'extern "C" '
    if ret == "complex<float>":
        arguments.append("float*")
        arguments.append("float*")
        ret = "void"
    if ret == "Vec" or ret == "vector<float>":
        arguments.append("float*")
        arguments.append("int*")
        ret = "void"
    if ret == "vector<Vec>" or ret == "vector<vector<float>>":
        arguments.append("float*")
        arguments.append("int*")
        arguments.append("int*")
        ret = "void"

    a = arg_cpp_to_c(ret, name)
    line += str(a) + " "
    line += name + "_export" + str(num)
    line += "("
    i = 1
    name = name.replace("_operator", "")
    while i < len(arguments):
        a = arg_cpp_to_c(arguments[i], name)
        if len(a) == 2:
            del arguments[i]
            arguments.insert(i, a[0])
            arguments.insert(i + 1, a[1])
            a = a[0]
        line += str(a) + " " + letters[i - 1]
        if i != len(arguments) - 1:
            line += ", "
        i += 1
    line += ") {"
    return line


def create_func(name, args, ret, num):
    lines = []
    i = 1
    if "operator" in name:
        if "complex" not in args[0]:
            line = "    return a->operator()("
        else:
            line = "    complex<float> r = a->operator()("
        i += 1
    else:
        if "ptr" in args[0]:
            line = "    return new " + name + "("
        else:
            line = "    return " + name + "("
    ind = i - 1
    while i < len(args):
        newline = type_convert_lines(args[i], i)
        lines.extend(newline)
        if len(newline) == 0:
            line += letters[ind]
        else:
            line += variables[args[i]]
            ind += 1
        if i != len(args) - 1:
            line += ", "
        ind += 1
        i += 1
    line += ");\n"
    if "complex" in args[0]:
        i = len(args) - 1
        if "Vec" in args:
            i += 1
        line += "    *" + letters[i] + " = real(r);\n"
        line += "    *" + letters[i + 1] + " = imag(r);\n"
    line += "}"
    lines.append(line)
    return lines


def write_export_include(INPUTS, lines):
    include_list = []
    for x in INPUTS:
        include_list.append('#include "../../' + x + '"\n')
    lines = base.remove_lines_between_phrases(
        lines, "// Begin include", "// End include"
    )
    lines = base.add_lines_between_phrases(
        lines, include_list, "// Begin include", "// End include"
    )
    return lines


""" Class Creation -------------
    elif "_" in name:
        var = str(name.split("_", 2)[-1])
        if "vector<vector<float>>" in args[0]:
            line = (
                "    *d = a->"
                + var
                + ".size();\n    for (int i = 0; i < *d; i++) {\n        c[i] = a->"
                + var
                + "[i].size();\n        for (int j = 0; j < c[i]; j++) {\n            b[i * c[i] + j] = a->"
                + var
                + "[i][j];\n        }\n    }\n"
            )
        elif "vector<Vec>" in args[0]:
            line = (
                "    *d = a->"
                + var
                + ".size();\n    for (int i = 0; i < *d; i++) {\n        c[i] = a->"
                + var
                + "[i].dimension;\n        for (int j = 0; j < c[i]; j++) {\n            b[i * c[i] + j] = a->"
                + var
                + "[i](j);\n        }\n    }\n"
            )
        elif "vector<float>" in args[0]:
            line = (
                "    *d = a->"
                + var
                + ".size();\n    for (int i = 0; i < *d; i++) {\n        b[i] = a->"
                + var
                + "[i];\n        }\n"
            )
        elif "Vec" in args[0]:
            line = (
                "    for (int i = 0; i < a->"
                + var
                + ".dimension; i++) {\n        b[i] = a->"
                + var
                + "(i);\n        c[i] = a->"
                + var
                + ".dimension;\n    }\n"
            )
        else:
            line = "    return a->" + var + ";\n"
"""


def write_export_func(INPUTS, lines):
    newlines = []
    for funcs in INPUTS.items():
        seen = defaultdict(int)
        i = 0
        name = funcs[0]
        args = funcs[1]["args"]
        ret = funcs[1]["return"]
        i = seen[name]
        seen[name] += 1
        line = create_declaration(name, args, ret, i) + "\n"
        func_lines = create_func(name, args, ret, i)
        for l in func_lines:
            line += l + "\n"
        prev_name = name
        newlines.append(line)
    lines = base.remove_lines_between_phrases(
        lines, "// Begin functions", "// End functions"
    )
    lines = base.add_lines_between_phrases(
        lines, newlines, "// Begin functions", "// End functions"
    )
    return lines


def write_export_funcs(INPUTS, FILES):
    input = copy.copy(INPUTS)
    with open("exports/cpp_export.cpp", "r") as file:
        lines = file.readlines()
    lines = write_export_include(FILES, lines)
    lines = write_export_func(input, lines)
    # with open("exports/cpp_export.cpp", "w") as file:
    #    file.writelines(lines)


def write_export_classes(INPUTS):

