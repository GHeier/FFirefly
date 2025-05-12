import re
from collections import defaultdict


def parse_args(arg_str):
    args = []
    if not arg_str.strip():
        return args

    # Handle comma-separated arguments with optional defaults
    for arg in re.split(r",(?![^()]*\))", arg_str):  # handle possible templates
        arg = arg.strip()
        if not arg:
            continue

        # Remove default values
        arg = arg.split("=")[0].strip()

        # Remove variable name, just keep type
        parts = arg.rsplit(" ", 1)
        if len(parts) == 2:
            arg_type = parts[0].strip()
        else:
            arg_type = parts[0].strip()  # handle pointer/ref cases
        args.append(arg_type)
    return args


def parse_functions(content):
    function_dict = defaultdict(list)
    function_pattern = re.compile(
        r"([a-zA-Z_][\w:<>\s*&]+?)\s+(operator[^\s(]*|[a-zA-Z_][\w]*)\s*\(([^)]*)\)\s*;"
    )

    for match in function_pattern.finditer(content):
        return_type = match.group(1).strip()
        name = match.group(2).strip()
        args = parse_args(match.group(3).strip())
        function_dict[name].append({"args": args, "return": return_type})

    return function_dict


def parse_classes(content):
    class_pattern = re.compile(r"class\s+(\w+)\s*{([^}]*)};", re.DOTALL)
    class_vars = {}
    class_methods = {}

    for class_match in class_pattern.finditer(content):
        class_name = class_match.group(1)
        body = class_match.group(2)

        vars = {}
        methods = defaultdict(list)

        for line in body.splitlines():
            line = line.strip().strip(";")
            if not line or line.startswith("//"):
                continue

            # Match constructor
            ctor_match = re.match(rf"({class_name})\s*\(([^)]*)\)", line)
            if ctor_match:
                name = ctor_match.group(1).strip()
                args = parse_args(ctor_match.group(2).strip())
                methods[name].append({"args": args, "return": ""})
                continue

            # Match operator() functions
            op_match = "operator" in line
            if op_match:
                return_type = line.split(" ")[0]
                args = parse_args(line.split("(")[-1][:-1])
                methods["operator"].append({"args": args, "return": return_type})
                continue

            # Match regular methods
            func_match = re.match(
                r"([a-zA-Z_][\w:<>\s*&]+?)\s+([a-zA-Z_][\w]*)\s*\(([^)]*)\)", line
            )
            if func_match:
                return_type = func_match.group(1).strip()
                name = func_match.group(2).strip()
                args = parse_args(func_match.group(3).strip())
                methods[name].append({"args": args, "return": return_type})
            else:
                # Skip lines with parentheses — they're function/method declarations
                if "(" in line or ")" in line:
                    continue

                # Match multi-variable lines like "float x, y, z;" or "int a = 1, b = 2;"
                var_decl_match = re.match(
                    r"([a-zA-Z_:<>0-9\s*&]+?)\s+([a-zA-Z_][\w\s,=]*)", line
                )
                if var_decl_match:
                    var_type = var_decl_match.group(1).strip()
                    var_names = [
                        v.strip().split("=")[0].strip()
                        for v in var_decl_match.group(2).split(",")
                    ]
                    for var in var_names:
                        if var:
                            vars[var] = var_type

        class_vars[class_name] = vars
        class_methods[class_name] = methods

    return class_vars, class_methods


def parse_hpp_file(filename):
    with open(filename, "r") as f:
        content = f.read()

    # Remove comments
    content = re.sub(r"//.*?$|/\*.*?\*/", "", content, flags=re.DOTALL | re.MULTILINE)

    class_vars, class_methods = parse_classes(content)
    global_funcs = parse_functions(content)

    # Flatten method names with index for overloads
    flattened_funcs = {}
    for name, overloads in global_funcs.items():
        if len(overloads) == 1:
            flattened_funcs[name] = overloads[0]
        else:
            for i, func in enumerate(overloads, 1):
                flattened_funcs[f"{name}{i}"] = func

    return flattened_funcs, class_vars, class_methods


def mark_constructors_as_ptr(class_methods):
    for class_name, methods in class_methods.items():
        if class_name in methods:
            for m in methods[class_name]:
                m["return"] = "ptr"


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def class_methods_to_funcs(class_methods):
    funcs = {}
    for cname, methods in class_methods.items():
        for name, overloads in methods.items():
            for i, func in enumerate(overloads, 1):
                label = f"{cname}_{name}"
                if name == cname:
                    label = f"{name}"
                funcs[label] = func
    return funcs


def global_funcs_remove_doubles(global_funcs, classes, class_methods):
    """
    Removes global functions that are identical (name, args, return) to class methods.
    """
    filtered_funcs = {}

    for func_name, func_data in global_funcs.items():
        # Normalize to list of overloads
        if isinstance(func_data, dict):
            func_data = [func_data]

        is_duplicate = False
        for class_name, methods in class_methods.items():
            if func_name in methods:
                class_overloads = methods[func_name]
                for f in func_data:
                    for m in class_overloads:
                        if f["args"] == m["args"] and f["return"] == m["return"]:
                            is_duplicate = True
                            break
                    if is_duplicate:
                        break
            if is_duplicate:
                break

        if not is_duplicate:
            # Store using consistent list format
            filtered_funcs[func_name] = func_data

    return filtered_funcs


def combine_function_dicts(global_funcs, class_methods, prefix_class=False):
    """
    Combine global functions and class methods into a single function dictionary.

    Args:
        global_funcs: dict of global function names -> list of overloads
        class_methods: dict of class names -> dict of method names -> list of overloads
        prefix_class: if True, prefix class methods with ClassName::MethodName to avoid conflicts

    Returns:
        Dict[str, List[Dict[str, Any]]]
    """
    combined = {}

    # Add global functions
    for name, overloads in global_funcs.items():
        combined[name] = overloads

    # Add class methods
    for class_name, methods in class_methods.items():
        for method_name, overloads in methods.items():
            key = f"{class_name}::{method_name}" if prefix_class else method_name
            if key in combined:
                combined[key].extend(overloads)
            else:
                combined[key] = overloads.copy()

    return combined


def print_all(global_funcs, classes, class_methods):
    print("Global Functions\n=======================\n")
    for func, args in global_funcs.items():
        print(func)
        print("    arg: ", args["args"])
        print("    ret: ", args["return"])
    print("Classes\n=======================\n")
    for cname, cvars in classes.items():
        print(cname)
        for vname, vtype in cvars.items():
            print(vname, ": ", vtype)
    print("Class Methods\n=======================\n")
    for name, funcs in class_methods.items():
        print(name)
        for func, args in funcs.items():
            for a in args:
                print(f"   {a['return']} {func}({a["args"]})")


from write import write_export, write_import
from write.base import add_lines_between_phrases, remove_lines_between_phrases

FILES = {
    "objects/vec.hpp",
    # "hamiltonian/band_structure.hpp",
    # "objects/CMField/bands.hpp",
}

#include_lines = write_export.write_export_include(FILES)
global_func_lines = []
class_func_lines = []
class_get_lines = []
for file in FILES:
    filename = "../" + file
    global_funcs, classes, class_methods = parse_hpp_file(filename)
    print_all(global_funcs, classes, class_methods)
    #global_funcs = global_funcs_remove_doubles(global_funcs, classes, class_methods)
    #global_funcs = combine_function_dicts(global_funcs, class_methods)
    #mark_constructors_as_ptr(class_methods)
    #class_funcs = class_methods_to_funcs(class_methods)

    #write_import.write_import(global_funcs, classes)

    #global_func_lines += write_export.write_export_class_func(
    #    global_funcs, list(class_methods)[0]
    #)
    ## class_func_lines += write_export.write_export_class_func(
    ##    class_funcs, list(class_methods)[0]
    ## )
    #class_get_lines += write_export.write_class_var_gets(classes)

#with open("exports/cpp_export.cpp", "r") as f:
#    lines = f.readlines()
#
#
## Include
#lines = remove_lines_between_phrases(lines, "// Begin include", "// End include")
#lines = add_lines_between_phrases(
#    lines, include_lines, "// Begin include", "// End include"
#)
#
## Class Functions
#lines = remove_lines_between_phrases(
#    lines, "// Begin Class functions", "// End Class functions"
#)
#lines = add_lines_between_phrases(
#    lines, class_func_lines, "// Begin Class functions", "// End Class functions"
#)
#
## Class Variable Functions
#lines = remove_lines_between_phrases(lines, "// Begin Class gets", "// End Class gets")
#lines = add_lines_between_phrases(
#    lines, class_get_lines, "// Begin Class gets", "// End Class gets"
#)
#
## Global Functions
#lines = remove_lines_between_phrases(lines, "// Begin functions", "// End functions")
#lines = add_lines_between_phrases(
#    lines, global_func_lines, "// Begin functions", "// End functions"
#)
#
#
#with open("exports/cpp_export.cpp", "w") as f:
#    f.writelines(lines)
