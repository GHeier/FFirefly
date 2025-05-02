import ctypes
import numpy as np
from typing import Any


def create_python_function(INPUTS):
    lines = []

    type_map = {
        "int": "ctypes.c_int",
        "float": "ctypes.c_float",
        "double": "ctypes.c_double",
        "string": "ctypes.c_char_p",
        "ptr": "ctypes.c_void_p",
        "complex<float>": "ctypes.c_float * 2",
        "None": "None",
    }

    friendly_names = {
        "int": "int",
        "float": "float",
        "double": "float",
        "string": "str",
        "ptr": "Any",
        "Vec": "list[float]",
        "complex<float>": "complex",
    }

    return_type_map = {
        "None": "None",
        "float": "float",
        "double": "float",
        "int": "int",
    }

    classes = {}
    functions = []

    for file, func_list in INPUTS.items():
        for base_name, sig in func_list:
            ret_type, *arg_types = sig
            entry = {
                "name": base_name,
                "ret_type": ret_type,
                "arg_types": arg_types,
            }
            if ret_type == "ptr":
                if base_name not in classes:
                    classes[base_name] = {"constructors": [], "methods": []}
                classes[base_name]["constructors"].append(entry)
            elif "_operator" in base_name:
                class_name = base_name.split("_operator")[0]
                if class_name in classes:
                    classes[class_name]["methods"].append(entry)
            else:
                functions.append(entry)

    # Free functions
    for entry in functions:
        base_name = entry["name"]
        arg_types = entry["arg_types"]
        ret_type = entry["ret_type"]
        func_name = f"{base_name}_export0"

        expanded_arg_ctypes = []
        annotations = []
        call_args = []

        param_index = 0
        for t in arg_types:
            name = f"arg{param_index}"
            if t == "Vec":
                expanded_arg_ctypes.extend(
                    ["ctypes.POINTER(ctypes.c_float)", "ctypes.c_int"]
                )
                annotations.append(f"{name}: list[float]")
                call_args.append(
                    f"{name}_array.ctypes.data_as(ctypes.POINTER(ctypes.c_float))"
                )
                call_args.append(f"len({name})")
            else:
                expanded_arg_ctypes.append(type_map[t])
                annotations.append(f"{name}: {friendly_names[t]}")
                call_args.append(name)
            param_index += 1

        lines.append(f"lib.{func_name}.argtypes = [{', '.join(expanded_arg_ctypes)}]\n")
        lines.append(f"lib.{func_name}.restype = {type_map[ret_type]}\n")
        lines.append("\n")

        lines.append(
            f"def {base_name}({', '.join(annotations)}) -> {return_type_map.get(ret_type, 'Any')}:\n"
        )
        for i, t in enumerate(arg_types):
            if t == "Vec":
                lines.append(f"    arg{i}_array = np.array(arg{i}, dtype=np.float32)\n")
        if ret_type != "None":
            lines.append(f"    result = lib.{func_name}({', '.join(call_args)})\n")
            lines.append("    return result\n")
        else:
            lines.append(f"    lib.{func_name}({', '.join(call_args)})\n")
        lines.append("\n")

    # Class-based objects
    for class_name, defs in classes.items():
        constructors = defs["constructors"]
        methods = defs["methods"]

        for i, ctor in enumerate(constructors):
            func = f"{class_name}_export{i}"
            arg_ctypes = [type_map[t] for t in ctor["arg_types"]]
            lines.append(f"lib.{func}.argtypes = [{', '.join(arg_ctypes)}]\n")
            lines.append(f"lib.{func}.restype = ctypes.c_void_p\n")
        if constructors:
            lines.append("\n")

        for i, method in enumerate(methods):
            func = f"{class_name}_operator_export{i}"
            expanded = []
            for t in method["arg_types"]:
                if t == "Vec":
                    expanded.extend(["ctypes.POINTER(ctypes.c_float)", "ctypes.c_int"])
                else:
                    expanded.append(type_map[t])
            lines.append(
                f"lib.{func}.argtypes = [ctypes.c_void_p, {', '.join(expanded)}]\n"
            )
            lines.append(f"lib.{func}.restype = ctypes.c_float\n")
        if methods:
            lines.append("\n")

        lines.append(f"class {class_name}:\n")
        lines.append("    def __init__(self, filename=None):\n")
        if any(len(c["arg_types"]) == 0 for c in constructors):
            lines.append(f"        if filename is None:\n")
            lines.append(f"            self.ptr = lib.{class_name}_export0()\n")
        if any(
            len(c["arg_types"]) == 1 and c["arg_types"][0] == "string"
            for c in constructors
        ):
            lines.append(f"        else:\n")
            lines.append(
                f"            self.ptr = lib.{class_name}_export1(ctypes.c_char_p(filename.encode('utf-8')))\n"
            )
        lines.append("        if not self.ptr:\n")
        lines.append(
            f"            raise RuntimeError('Failed to initialize {class_name}')\n"
        )
        lines.append("\n")

        lines.append("    def __call__(self, *args):\n")
        for i, method in enumerate(methods):
            arg_types = method["arg_types"]
            cond = []
            pre = []
            call = []
            argc = 0

            for j, t in enumerate(arg_types):
                if t == "float":
                    cond.append(f"isinstance(args[{j}], (int, float))")
                    pre.append(f"        arg{argc} = ctypes.c_float(args[{j}])")
                    call.append(f"arg{argc}")
                    argc += 1
                elif t == "int":
                    cond.append(f"isinstance(args[{j}], int)")
                    pre.append(f"        arg{argc} = ctypes.c_int(args[{j}])")
                    call.append(f"arg{argc}")
                    argc += 1
                elif t == "Vec":
                    cond.append(f"isinstance(args[{j}], (list, tuple))")
                    pre.append(
                        f"        arg{argc} = (ctypes.c_float * len(args[{j}]))(*[float(x) for x in args[{j}]])"
                    )
                    pre.append(f"        arg{argc+1} = ctypes.c_int(len(args[{j}]))")
                    call.extend([f"arg{argc}", f"arg{argc+1}"])
                    argc += 2
                else:
                    cond.append("True")

            cond_str = " and ".join(cond)
            lines.append(f"        if len(args) == {len(arg_types)} and {cond_str}:\n")
            lines.extend(line + "\n" for line in pre)
            lines.append(
                f"            return lib.{class_name}_operator_export{i}(self.ptr, {', '.join(call)})\n"
            )
        lines.append("        raise TypeError('Invalid arguments to __call__')\n")
        lines.append("\n")

        lines.append("    def __del__(self):\n")
        lines.append("        try:\n")
        lines.append("            destroy = lib.destroy\n")
        lines.append("            destroy.argtypes = [ctypes.c_void_p]\n")
        lines.append("            destroy(self.ptr)\n")
        lines.append("        except AttributeError:\n")
        lines.append("            pass\n")
        lines.append("\n")

    return lines


# Python function wrapper generator
# def generate_python_function(name, arg_names):
#    args_str = ", ".join(arg_names)
#    call_args = ", ".join(arg_names)
#    return f"""
# def {name}_py({args_str}):
#    return {name}({call_args})
# """
#
#
# def create_function(INPUTS):
#    # Convert to flat list of functions
#    functions = []
#    for file, funcs in INPUTS.items():
#        for name, args in funcs:
#            functions.append(
#                {
#                    "name": f"{name}_export0",  # Assuming suffix _export0, adjust if needed
#                    "args": args,
#                    "return": (
#                        "float" if name == "epsilon" else "None"
#                    ),  # Basic heuristic
#                }
#            )
#
#    # Generate ctypes wrappers
#    for f in functions:
#        func = getattr(lib, f["name"])
#        func.argtypes = [type_map[arg] for arg in f["args"]]
#        func.restype = type_map[f["return"]]
#        globals()[f["name"]] = func
#
#    # Print example Python function wrappers
#    for f in functions:
#        arg_names = [f"arg{i}" for i in range(len(f["args"]))]
#        print(generate_python_function(f["name"], arg_names))
