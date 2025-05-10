import ctypes
import numpy as np
from typing import Any


def create_python_function(INPUTS):
    lines = []

    # === Type Maps ===
    type_map = {
        "int": "ctypes.c_int",
        "float": "ctypes.c_float",
        "double": "ctypes.c_double",
        "string": "ctypes.c_char_p",
        "ptr": "ctypes.c_void_p",
        "complex<float>": "ctypes.c_float * 2",  # special handling
        "None": "None",
    }

    friendly_names = {
        "int": "int",
        "float": "float",
        "double": "float",
        "string": "str",
        "ptr": "Any",
        "Vec": "list[float]",
        "vector<float>": "list[float]",
        "complex<float>": "complex",
        "vector<Vec>": "list[list[float]]",
        "vector<vector<float>>": "list[list[float]]",
    }

    return_type_map = {
        "None": "None",
        "float": "float",
        "double": "float",
        "int": "int",
        "complex<float>": "complex",
    }

    def parse_arg(arg):
        if "=" in arg:
            t, default = arg.split("=")
            return t.strip(), default.strip()
        return arg.strip(), None

    # Makes inputs into classes and functions
    # Classes are based on returning ptr, functions are added to them if they contain "_operator"
    def split_inputs(INPUTS):
        classes = {}
        functions = []

        for file, func_list in INPUTS.items():
            for base_name, sig in func_list:
                ret_type, *raw_args = sig
                parsed_args = [parse_arg(a) for a in raw_args]

                if "_var_" in base_name:
                    # e.g. Surface_var_faces -> class: Surface, var: faces
                    class_name, var_name = base_name.split("_var_", 1)
                    if class_name not in classes:
                        classes[class_name] = {
                            "constructors": [],
                            "methods": [],
                            "variables": [],
                        }
                    classes[class_name].setdefault("variables", []).append(
                        (var_name, ret_type)
                    )
                    continue

                entry = {
                    "name": base_name,
                    "ret_type": ret_type,
                    "arg_types": parsed_args,
                }

                if ret_type == "ptr":
                    if base_name not in classes:
                        classes[base_name] = {
                            "constructors": [],
                            "methods": [],
                            "variables": [],
                        }
                    classes[base_name]["constructors"].append(entry)
                elif "_operator" in base_name:
                    class_name = base_name.split("_operator")[0]
                    if class_name in classes:
                        classes[class_name]["methods"].append(entry)
                else:
                    functions.append(entry)

        return functions, classes

    def emit_free_function(entry):
        base_name = entry["name"]
        arg_types = entry["arg_types"]
        ret_type = entry["ret_type"]
        func_name = f"{base_name}_export0"

        expanded_arg_ctypes = []
        annotations = []
        call_args = []

        for i, (t, default) in enumerate(arg_types):
            name = f"arg{i}"
            if t == "Vec" or t == "vector<float>":
                expanded_arg_ctypes += [
                    "ctypes.POINTER(ctypes.c_float)",
                    "ctypes.c_int",
                ]
                annotations.append(f"{name}: list[float]")
                call_args += [
                    f"{name}_array.ctypes.data_as(ctypes.POINTER(ctypes.c_float))",
                    f"len({name})",
                ]
            elif t == "vector<Vec>" or t == "vector<vector<float>>":
                expanded_arg_ctypes += [
                    "ctypes.POINTER(ctypes.c_float)",
                    "ctypes.POINTER(ctypes.c_int)",
                    "ctypes.c_int",
                ]
                annotations.append(f"{name}: list[list[float]]")
                call_args.append(f"{name}_flat_arr")
                call_args.append(f"{name}_len_arr")
                call_args.append(f"len(len_{name})")
                lines.append(f"    flat_{name} = []\n")
                lines.append(f"    len_{name} = []\n")
                lines.append(f"    for vec in {name}:\n")
                lines.append(f"        flat_{name}.extend(vec)\n")
                lines.append(f"        len_{name}.append(len(vec))\n")
                lines.append(
                    f"    {name}_flat_arr = (ctypes.c_float * len(flat_{name}))(*flat_{name})\n"
                )
                lines.append(
                    f"    {name}_len_arr = (ctypes.c_int * len(len_{name}))(*len_{name})\n"
                )
            else:
                expanded_arg_ctypes.append(type_map[t])
                arg_decl = f"{name}: {friendly_names[t]}"
                if default:
                    arg_decl += f" = {default}"
                annotations.append(arg_decl)
                call_args.append(name)

        lines.append(f"lib.{func_name}.argtypes = [{', '.join(expanded_arg_ctypes)}]\n")
        lines.append(f"lib.{func_name}.restype = {type_map[ret_type]}\n\n")

        lines.append(
            f"def {base_name}({', '.join(annotations)}) -> {return_type_map.get(ret_type, 'Any')}:\n"
        )
        for i, (t, _) in enumerate(arg_types):
            if t == "Vec" or t == "vector<float>":
                lines.append(f"    arg{i}_array = np.array(arg{i}, dtype=np.float32)\n")
        if ret_type != "None":
            lines.append(f"    result = lib.{func_name}({', '.join(call_args)})\n")
            lines.append("    return result\n")
        else:
            lines.append(f"    lib.{func_name}({', '.join(call_args)})\n")
        lines.append("\n")

    def emit_constructor_bindings(class_name, constructors):
        for i, ctor in enumerate(constructors):
            func = f"{class_name}_export{i}"
            arg_ctypes = [type_map[t] for t, _ in ctor["arg_types"]]
            lines.append(f"lib.{func}.argtypes = [{', '.join(arg_ctypes)}]\n")
            lines.append("lib.func.restype = ctypes.c_void_p\n".replace("func", func))
        if constructors:
            lines.append("\n")

    def emit_method_bindings(class_name, methods):
        for i, method in enumerate(methods):
            func = f"{class_name}_operator_export{i}"
            ret_type = method["ret_type"]
            arg_types = method["arg_types"]

            if arg_types and arg_types[0][0] == "ptr":
                arg_types = arg_types[1:]

            arg_list = []

            for t, _ in arg_types:
                if t == "Vec" or t == "vector<float>":
                    arg_list += ["ctypes.POINTER(ctypes.c_float)", "ctypes.c_int"]
                else:
                    arg_list.append(type_map[t])

            if ret_type == "complex<float>":
                arg_list += [
                    "ctypes.POINTER(ctypes.c_float)",
                    "ctypes.POINTER(ctypes.c_float)",
                ]
            lines.append(
                f"lib.{func}.argtypes = [ctypes.c_void_p, {', '.join(arg_list)}]\n"
            )
            lines.append(
                f"lib.{func}.restype = {'None' if ret_type == 'complex<float>' else type_map[ret_type]}\n"
            )
        if methods:
            lines.append("\n")

    def emit_class(class_name, constructors, methods):
        lines.append(f"class {class_name}:\n")
        # Emit class variables from Surface_var_xxx
        for var_name, var_type in defs.get("variables", []):
            if var_type == "vector<Vec>" or var_type == "vector<vector<float>>":
                py_type = "list[list[float]]"
                init_val = "[]"
            elif var_type == "Vec" or var_type == "vector<float>":
                py_type = "list[float]"
                init_val = "[]"
            elif var_type in friendly_names:
                py_type = friendly_names[var_type]
                init_val = "None"
            else:
                py_type = "Any"
                init_val = "None"
            lines.append(f"    {var_name}: {py_type} = {init_val}\n")

        lines.append("    def __init__(self, filename=None):\n")
        if any(len(c["arg_types"]) == 0 for c in constructors):
            lines.append("        if filename is None:\n")
            lines.append(f"            self.ptr = lib.{class_name}_export0()\n")
        if any(
            len(c["arg_types"]) == 1 and c["arg_types"][0][0] == "string"
            for c in constructors
        ):
            lines.append("        else:\n")
            lines.append(
                f"            self.ptr = lib.{class_name}_export1(ctypes.c_char_p(filename.encode('utf-8')))\n"
            )

        for var_name, var_type in defs.get("variables", []):
            export_func = f"{class_name}_var_{var_name}_export0"
            lines.append(f"        # Load variable '{var_name}' from C++\n")
            lines.append(f"        count = ctypes.c_int()\n")
            lines.append(f"        lib.{export_func}.restype = None\n")
            lines.append(
                f"        lib.{export_func}.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_int), ctypes.c_int]\n"
            )
            lines.append(f"        # First call to get number of vectors\n")
            lines.append(
                f"        lib.{export_func}(self.ptr, None, None, ctypes.byref(count))\n"
            )
            lines.append(f"        n = count.value\n")
            lines.append(f"        # Allocate buffers\n")
            lines.append(f"        lens = (ctypes.c_int * n)()\n")
            lines.append(f"        total_len = 0\n")
            lines.append(
                f"        lib.{export_func}(self.ptr, None, lens, ctypes.c_int(n))\n"
            )
            lines.append(f"        for i in range(n):\n")
            lines.append(f"            total_len += lens[i]\n")
            lines.append(f"        buf = (ctypes.c_float * total_len)()\n")
            lines.append(
                f"        lib.{export_func}(self.ptr, buf, lens, ctypes.c_int(n))\n"
            )
            lines.append(f"        offset = 0\n")
            lines.append(f"        result = []\n")
            lines.append(f"        for i in range(n):\n")
            lines.append(
                f"            sub = [buf[offset + j] for j in range(lens[i])]\n"
            )
            lines.append(f"            result.append(sub)\n")
            lines.append(f"            offset += lens[i]\n")
            lines.append(f"        self.{var_name} = result\n\n")
        lines.append("        if not self.ptr:\n")
        lines.append(
            f"            raise RuntimeError('Failed to initialize {class_name}')\n\n"
        )

        lines.append("    def __call__(self, *args):\n")
        for i, method in enumerate(methods):
            ret_type = method["ret_type"]
            arg_types = method["arg_types"]
            if arg_types and arg_types[0][0] == "ptr":
                arg_types = arg_types[1:]

            # Count required arguments
            required = 0
            for t, default in reversed(arg_types):
                if default is None:
                    break
                required += 1
            min_args = len(arg_types) - required
            max_args = len(arg_types)

            cond = [f"len(args) >= {min_args}", f"len(args) <= {max_args}"]
            pre, call = [], []
            argc = 0

            for j, (t, default) in enumerate(arg_types):
                name = f"arg{argc}"
                if t == "float":
                    if default:
                        pre.append(
                            f"{name} = ctypes.c_float(args[{j}]) if len(args) > {j} else ctypes.c_float({default})"
                        )
                    else:
                        cond.append(f"isinstance(args[{j}], (int, float))")
                        pre.append(f"{name} = ctypes.c_float(args[{j}])")
                    call.append(name)
                    argc += 1
                elif t == "int":
                    if default:
                        pre.append(
                            f"{name} = ctypes.c_int(args[{j}]) if len(args) > {j} else ctypes.c_int({default})"
                        )
                    else:
                        cond.append(f"isinstance(args[{j}], int)")
                        pre.append(f"{name} = ctypes.c_int(args[{j}])")
                    call.append(name)
                    argc += 1
                elif t == "string":
                    if default:
                        pre.append(
                            f"{name} = ctypes.c_char_p(args[{j}].encode('utf-8')) if len(args) > {j} else ctypes.c_char_p(b{default})"
                        )
                    else:
                        cond.append(f"isinstance(args[{j}], str)")
                        pre.append(
                            f"{name} = ctypes.c_char_p(args[{j}].encode('utf-8'))"
                        )
                    call.append(name)
                    argc += 1
                elif t == "Vec" or t == "vector<float>":
                    if default:
                        cond.append(f"isinstance(args[{j}], (list, tuple))")
                        pre.append(
                            f"{name} = (ctypes.c_float * len(args[{j}]))(*[float(x) for x in args[{j}]]) if len(args) > {j} else ctypes.c_float * len(args[{j}]))(*{default})"
                        )
                    else:
                        pre.append(
                            f"{name} = (ctypes.c_float * len(args[{j}]))(*[float(x) for x in args[{j}]])"
                        )
                    pre.append(f"{name}_len = ctypes.c_int(len(args[{j}]))")
                    call.extend([name, f"{name}_len"])
                    argc += 2

            lines.append(
                f"        # Overload for args={len(arg_types)}, required={min_args}\n"
            )
            lines.append(f"        if {' and '.join(cond)}:\n")
            if ret_type == "complex<float>":
                lines.append("            real = ctypes.c_float()\n")
                lines.append("            imag = ctypes.c_float()\n")
                lines.extend([f"            {line}\n" for line in pre])
                lines.append(
                    f"            lib.{class_name}_operator_export{i}(self.ptr, {', '.join(call)}, ctypes.byref(real), ctypes.byref(imag))\n"
                )
                lines.append("            return complex(real.value, imag.value)\n")
            else:
                lines.extend([f"            {line}\n" for line in pre])
                lines.append(
                    f"            return lib.{class_name}_operator_export{i}(self.ptr, {', '.join(call)})\n"
                )
        lines.append("        raise TypeError('Invalid arguments to __call__')\n\n")

        lines.append("    def __del__(self):\n")
        lines.append("        try:\n")
        lines.append("            destroy = lib.destroy\n")
        lines.append("            destroy.argtypes = [ctypes.c_void_p]\n")
        lines.append("            destroy(self.ptr)\n")
        lines.append("        except AttributeError:\n")
        lines.append("            pass\n\n")

    # === Main Generation ===
    functions, classes = split_inputs(INPUTS)

    for f in functions:
        emit_free_function(f)

    for class_name, defs in classes.items():
        emit_constructor_bindings(class_name, defs["constructors"])
        emit_method_bindings(class_name, defs["methods"])
        emit_class(class_name, defs["constructors"], defs["methods"])

    return lines
