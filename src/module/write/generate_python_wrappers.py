import ctypes
import numpy as np
from typing import Any

letters = "abcdefghijklmnopq"


def create_python_function(global_funcs, classes):
    lines = []

    # === Type Maps ===
    type_map = {
        "bool": "ctypes.c_bool",
        "int": "ctypes.c_int",
        "float": "ctypes.c_float",
        "double": "ctypes.c_double",
        "string": "ctypes.c_char_p",
        "ptr": "ctypes.c_void_p",
        "Vec": "ctypes.c_void_p",
        "const Vec": "ctypes.c_void_p",
        "vector<float>": "ctypes.c_void_p",
        "const float": "ctypes.c_float",
        "complex<Vec>": "ctypes.c_void_p * 2",
        "const complex<Vec>": "ctypes.c_void_p * 2",
        "complex<float>": "ctypes.c_float * 2",  # special handling
        "None": "None",
    }

    friendly_names = {
        "bool": "bool",
        "int": "int",
        "float": "float",
        "const float": "float",
        "double": "float",
        "string": "str",
        "ptr": "any",
        "Vec": "Vec",
        "const Vec": "Vec",
        "vector<float>": "list[float]",
        "complex<float>": "complex",
        "complex<Vec>": "(Vec, Vec)",
        "const complex<Vec>": "(Vec, Vec)",
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

    def split_inputs(global_funcs, class_fields):
        classes = {}
        functions = []

        for name, overloads in global_funcs.items():
            for overload in overloads:
                if overload["return"] == "explicit":
                    overload["return"] = "ptr"
                ret_type = overload["return"]
                arg_types = [(arg.strip(), None) for arg in overload["args"]]

                entry = {
                    "name": name,
                    "ret_type": ret_type,
                    "arg_types": arg_types,
                }

                if name in class_fields:
                    # This is either a constructor or method
                    if ret_type == "ptr":
                        if name not in classes:
                            classes[name] = {
                                "constructors": [],
                                "methods": [],
                                "fields": class_fields[name],
                            }
                        classes[name]["constructors"].append(entry)
                    else:
                        # This might be a static method returning the class
                        functions.append(entry)
                elif "_operator" in name:
                    class_name = name.split("_operator")[0]
                    if class_name in classes:
                        classes[class_name]["methods"].append(entry)
                    else:
                        classes[class_name] = {
                            "constructors": [],
                            "methods": [entry],
                            "fields": {},
                        }
                else:
                    functions.append(entry)

        # Include class field-only definitions (no constructor yet)
        for cname, fields in class_fields.items():
            if cname not in classes:
                classes[cname] = {"constructors": [], "methods": [], "fields": fields}

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
            f"def {base_name}({', '.join(annotations)}) -> {return_type_map.get(ret_type, 'any')}:\n"
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

    def emit_class(class_name, constructors, methods, fields):
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
                py_type = "any"
                init_val = "None"
            lines.append(f"    {var_name}: {py_type} = {init_val}\n")

        arglines = []
        callines = []
        for i in range(len(constructors)):
            c = constructors[i]
            alength = len(c["arg_types"])

            next_constructor = False
            argline = "if len(args) == "
            calline = ""
            if i > 0:
                argline = "el" + argline
            argline += str(alength)
            if alength > 0:
                for ai in range(alength):
                    if c["arg_types"][ai][0] not in friendly_names:
                        next_constructor = True
                        break
                    argline += f" and isinstance(args[{ai}], {friendly_names[c['arg_types'][ai][0]]}):\n"
                    calline += letters[ai] + ", "
            else:
                argline += ":\n"
            if next_constructor:
                continue
            calline = calline[:-2]
            arglines.append(argline)
            callines.append(calline)

        field_lines = []
        for vname, vtype in fields.items():
            fline = (
                f"        self.{vname} = lib.{class_name}_{vname}_export0(self.ptr)\n"
            )
            field_lines.append(fline)

        lines.append(f"    def __init__(self, *args):\n")
        for i in range(len(arglines)):
            l = arglines[i]
            calline = callines[i]
            lines.append(f"        {l}")
            lines.append(
                f"            self.ptr = lib.{class_name}_export{i}({calline})\n"
            )
            lines.append("            if not self.ptr:\n")
            lines.append(
                f"                raise RuntimeError('Failed to initialize {class_name}')\n\n"
            )
        for f in field_lines:
            lines.append(f)

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
    functions, classes = split_inputs(global_funcs, classes)

    for f in functions:
        emit_free_function(f)

    for class_name, defs in classes.items():
        emit_constructor_bindings(class_name, defs["constructors"])
        emit_method_bindings(class_name, defs["methods"])
        emit_class(class_name, defs["constructors"], defs["methods"], defs["fields"])

    return lines
