from . import base
import copy
from collections import defaultdict

letters = "abcdefghijklmnopq"
variables = {"Vec": "v"}
special_types = {"Vec"}


def arg_cpp_to_c(arg, name):
    if arg == "Vec":
        return "Vec*"
    if arg == "vector":
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
            + letters[i]
            + ", "
            + letters[i + 1]
            + ");"
        )
    return lines


def filter_default_args(args):
    arguments = []
    for x in args:
        arguments.append(x.split("=")[0])
    return arguments


def create_class_declaration(classname, name, args, ret, num):
    arguments = filter_default_args(args)

    line = 'extern "C" '
    if ret == "complex<float>":
        arguments.append("float*")
        arguments.append("float*")
        ret = "void"
    if ret == classname:
        ret = classname + "*"
    if ret == "vector<float>":
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
    i = 0
    name = name.replace("_operator", "")
    while i < len(arguments):
        a = arg_cpp_to_c(arguments[i], name)
        if len(a) == 2:
            del arguments[i]
            arguments.insert(i, a[0])
            arguments.insert(i + 1, a[1])
            a = a[0]
        line += str(a) + " " + letters[i]
        if i != len(arguments) - 1:
            line += ", "
        i += 1
    line += ") {"
    return line


def create_declaration(name, args, ret, num):
    arguments = filter_default_args(args)

    line = 'extern "C" '
    if ret == "complex<float>":
        arguments.append("float*")
        arguments.append("float*")
        ret = "void"
    if ret == "Vec":
        ret = "Vec*"
    if ret == "vector<float>":
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
    i = 0
    name = name.replace("_operator", "")
    while i < len(arguments):
        a = arg_cpp_to_c(arguments[i], name)
        if len(a) == 2:
            del arguments[i]
            arguments.insert(i, a[0])
            arguments.insert(i + 1, a[1])
            a = a[0]
        line += str(a) + " " + letters[i]
        if i != len(arguments) - 1:
            line += ", "
        i += 1
    line += ") {"
    return line


def create_var_declaration(name, args, ret, num):
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
    i = 0
    name = name.replace("_operator", "")
    name = name.replace("__", "_")
    name, vname = name.split("_")
    while i < len(arguments):
        a = arg_cpp_to_c(arguments[i], name)
        if len(a) == 2:
            del arguments[i]
            arguments.insert(i, a[0])
            arguments.insert(i + 1, a[1])
            a = a[0]
        line += str(a) + " " + letters[i]
        if i != len(arguments) - 1:
            line += ", "
        i += 1
    line += ") {"
    return line


def create_class_func(classname, name, args, ret, num):
    lines = []
    i = 0
    i += len(args) > 0 and args[0] == classname + "*"
    if "operator" in name:
        if "complex" not in ret:
            line = "    return a->operator()("
        else:
            line = "    complex<float> r = a->operator()("
        i += 1
    else:
        if "ptr" in ret:
            line = "    return new " + name + "("
        elif ret in special_types:
            line = f"    {ret}* result = new {ret}({name}("
            if ret == name:
                line = f"    {ret}* result = new {ret}("
        elif "string" in ret:
            line = "    return strdup(" + name + "("
        elif "vector<float>" in ret:
            line = "    vector<float> result = " + name + "("
        else:
            line = "    return a->" + name.split("_")[-1] + "("
    if len(args) > 0 and args[0] == ret:
        line = f"    {ret}* result = new {ret}(a->{name.split("_")[-1]}("
    ind = i
    while i < len(args):
        # newline = type_convert_lines(args[i], i)
        # lines.extend(newline)
        # if len(newline) == 0:
        #    line += letters[ind]
        # else:
        #    line += variables[args[i]]
        #    ind += 1
        if args[i] == classname:
            line += "*"
        line += letters[ind]
        if i != len(args) - 1:
            line += ", "
        ind += 1
        i += 1
    if classname == ret and ret != name:
        line += ")"
    if ret == "string":
        line += ").c_str());\n"
    else:
        line += ");\n"
    if "complex" in ret:
        i = len(args) - 1
        if "Vec" in args:
            i += 1
        line += "    *" + letters[i] + " = real(r);\n"
        line += "    *" + letters[i + 1] + " = imag(r);\n"
    if classname == ret:
        line += f"    return result;\n"
    line += "}"
    lines.append(line)
    return lines


def create_func(name, args, ret, num):
    lines = []
    i = 0
    if "operator" in name:
        if "complex" not in ret:
            line = "    return a->operator()("
        else:
            line = "    complex<float> r = a->operator()("
        i += 1
    else:
        if "ptr" in ret:
            line = "    return new " + name + "("
        elif ret in special_types:
            line = f"    {ret}* result = new {ret}({name}("
        elif "string" in ret:
            line = "    return strdup(" + name + "("
        elif "vector<float>" in ret:
            line = "    vector<float> result = " + name + "("
        else:
            line = "    return " + name + "("
    ind = i
    while i < len(args):
        # newline = type_convert_lines(args[i], i)
        # lines.extend(newline)
        # if len(newline) == 0:
        #    line += letters[ind]
        # else:
        #    line += variables[args[i]]
        #    ind += 1
        if args[i] in special_types:
            line += "*"
        line += letters[ind]
        if i != len(args) - 1:
            line += ", "
        ind += 1
        i += 1
    if ret == "string":
        line += ").c_str());\n"
    else:
        if ret in special_types:
            line += ")"
        line += ");\n"
    if "complex" in ret:
        i = len(args) - 1
        if "Vec" in args:
            i += 1
        line += "    *" + letters[i] + " = real(r);\n"
        line += "    *" + letters[i + 1] + " = imag(r);\n"
    if ret == "vector<float>":
        line += f"    vector_to_ptr(result, {letters[i]}, {letters[i+1]});\n"
    if ret in special_types:
        line += f"    return result;\n"
    line += "}"
    lines.append(line)
    return lines


def write_export_include(INPUTS):
    include_list = []
    for x in INPUTS:
        include_list.append('#include "../../' + x + '"\n')
    return include_list


def write_export_func(INPUTS):
    newlines = []
    seen = defaultdict(int)
    for funcs in INPUTS.items():
        i = 0
        name = funcs[0]
        args = funcs[1][0]["args"]
        ret = funcs[1][0]["return"]
        if ret == "explicit":
            ret = ""
        i = seen[name]
        seen[name] += 1
        line = create_declaration(name, args, ret, i) + "\n"
        func_lines = create_func(name, args, ret, i)
        for l in func_lines:
            line += l + "\n"
        newlines.append(line)
    return newlines


def write_export_class_func(INPUTS, classname):
    newlines = []
    seen = defaultdict(int)
    i = 0
    for funcs in INPUTS.items():
        name = funcs[0]
        for f in funcs[1]:
            args = f["args"]
            ret = f["return"]
            if ret == "explicit":
                ret = "Vec"
            if name != classname and (len(args) == 0 or args[0] != classname):
                args.insert(0, classname + "*")
            i = seen[name]
            seen[name] += 1
            line = create_class_declaration(classname, name, args, ret, i) + "\n"
            func_lines = create_class_func(classname, name, args, ret, i)
            for l in func_lines:
                line += l + "\n"
            newlines.append(line)
    return newlines


def create_var_func(vname, args, ret):
    newlines = []
    newlines.append(f"    return a->{vname};\n" + "}")
    return newlines


def write_class_var_gets(INPUTS):
    newlines = []
    for funcs in INPUTS.items():
        name = funcs[0]
        for vname, vtype in funcs[1].items():
            nname = name + "_" + vname
            args = ["ptr"]
            ret = vtype
            line = create_var_declaration(nname, args, ret, 0) + "\n"
            func_lines = create_var_func(vname, args, ret)
            for l in func_lines:
                line += l + "\n"
            newlines.append(line)
    return newlines
