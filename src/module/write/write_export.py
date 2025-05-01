from . import base
import copy

letters = "abcdefghijklmnopqrstu"
variables = {"Vec": "v"}


def write_export_include(INPUTS, lines):
    include_list = []
    for x in INPUTS:
        include_list.append(x)
    remove_lines_between_phrases(lines, "// Begin include", "// End include")
    add_lines_between_phrases(lines, include_list, "// Begin include", "// End include")
    return lines


def arg_cpp_to_c(arg, name):
    if arg == "Vec":
        return ["float*", "int"]
    if arg == "vector":
        return [""]
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
    if a == "vector":
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


def create_declaration(name, args, num):
    arguments = copy.copy(args)
    line = 'extern "C" '
    a = arg_cpp_to_c(arguments[0], name)
    line += a + " "
    line += name + "_export" + str(num)
    line += "("
    i = 1
    while i < len(arguments):
        a = arg_cpp_to_c(arguments[i], name)
        if len(a) == 2:
            del arguments[i]
            arguments.insert(i, a[0])
            arguments.insert(i + 1, a[1])
            a = a[0]
        line += a + " " + letters[i - 1]
        if i != len(arguments) - 1:
            line += ", "
        i += 1
    line += ") {"
    return line


def create_func(name, args, num):
    lines = []

    line = "    return " + name + "("
    for i in range(1, len(args)):
        newline = type_convert_lines(args[i], i)
        lines.extend(newline)
        if len(newline) == 0:
            line += letters[i - 1]
        else:
            line += variables[args[i]]
        if i != len(args) - 1:
            line += ", "
    line += ");\n}"
    lines.append(line)
    return lines


def write_export_func(INPUTS, lines):
    for file, func in INPUTS.items():
        i = 0
        prev_name = ""
        for name, args in func.items():
            line = create_declaration(name, args, i) + "\n"
            func_lines = create_func(name, args, i)
            for l in func_lines:
                line += l + "\n"
            if prev_name == name:
                i += 1
            else:
                i = 0
            prev_name = name
            print(line)
            lines.append(line)
    return lines


def write_export(INPUTS):
    with open("exports/cpp_export.cpp", "r") as file:
        lines = file.readlines()
    write_export_func(INPUTS, lines)
