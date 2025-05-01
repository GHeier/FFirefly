from . import base

letters = "abcdefghijklmnopqrstu"


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


def type_convert_lines(args):
    lines = []
    for i in range(len(args)):
        a = args[i]
        if a == "Vec":
            lines.append("Vec v(" + letters[i] + ", " + letters[i + 1] + ");")
        if a == "vector":
            lines.append("Vec v(" + letters[i] + ", " + letters[i + 1] + ");")


def create_declaration(name, args, num):
    line = 'extern "C" '
    a = arg_cpp_to_c(args[0], name)
    line += a + " "
    line += name + "_export" + str(num)
    line += "("
    for i in range(1, len(args)):
        a = arg_cpp_to_c(args[i], name)
        if len(a) == 2:
            args.insert(i, a[1])
        line += a + " " + letters[i - 1]
        if i != len(args) - 1:
            line += ", "
    line += ") {"
    return line


def write_export_func(INPUTS, lines):
    for file, func in INPUTS.items():
        i = 0
        prev_name = ""
        for name, args in func.items():
            line = create_declaration(name, args, i)
            if prev_name == name:
                i += 1
            else:
                i = 0
            prev_name = name
            print(line)
            lines.append(line)
    return lines


def write_export(INPUTS):
    with open("../exports/cpp_export.cpp", "r") as file:
        lines = file.readlines()
    write_export_func(INPUTS, lines)
