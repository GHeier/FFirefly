def create_julia_function(INPUTS):
    lines = []
    exports = []

    type_map = {
        "int": "Cint",
        "float": "Float32",
        "double": "Float64",
        "string": "Cstring",
        "ptr": "Ptr{Cvoid}",
        "complex<float>": "ComplexF32",
        "None": "Nothing",
    }

    julia_types = {
        "int": "Int",
        "float": "Float32",
        "double": "Float64",
        "string": "String",
        "ptr": "Ptr{Cvoid}",
        "Vec": "Vector{Float64}",
        "complex<float>": "ComplexF32",
    }

    def parse_arg(arg):
        if "=" in arg:
            t, default = arg.split("=")
            return t.strip(), default.strip()
        return arg.strip(), None

    def split_inputs(INPUTS):
        classes = {}
        functions = []
        for file, func_list in INPUTS.items():
            for base_name, sig in func_list:
                ret_type, *raw_args = sig
                parsed_args = [parse_arg(a) for a in raw_args]
                entry = {
                    "name": base_name,
                    "ret_type": ret_type,
                    "arg_types": parsed_args,
                }
                if ret_type == "ptr":
                    if base_name not in classes:
                        classes[base_name] = {"constructors": [], "methods": []}
                    classes[base_name]["constructors"].append(entry)
                    exports.append(base_name)
                elif "_operator" in base_name:
                    class_name = base_name.split("_operator")[0]
                    if class_name in classes:
                        classes[class_name]["methods"].append(entry)
                else:
                    functions.append(entry)
                    exports.append(base_name)
        return functions, classes

    def emit_free_function(entry):
        name = entry["name"]
        ret = entry["ret_type"]
        args = entry["arg_types"]

        jl_args = []
        ccall_types = []
        ccall_vals = []
        pre = []

        for i, (t, default) in enumerate(args):
            var = f"arg{i}"
            if t == "Vec":
                jl_args.append(f"{var}::Vector{{Float64}}")
                pre.append(f"    new{var} = Float32.({var})")
                ccall_types.extend(["Ptr{Float32}", "Cint"])
                ccall_vals.extend([f"new{var}", f"length({var})"])
            else:
                jl_args.append(
                    f"{var}::{julia_types[t]}" + (f"={default}" if default else "")
                )
                ccall_types.append(type_map[t])
                ccall_vals.append(var)

        lines.append(f"function {name}({', '.join(jl_args)})\n")
        lines.extend(pre)
        lines.append(
            f"\n    return ccall((:{name}_export0, libfly), {type_map[ret]}, ({', '.join(ccall_types)}), {', '.join(ccall_vals)})\n"
        )
        lines.append("end\n\n")

    def emit_class(class_name, constructors, methods):
        lines.append(f"mutable struct {class_name}\n    ptr::Ptr{{Cvoid}}\nend\n\n")

        for i, ctor in enumerate(constructors):
            args = ctor["arg_types"]
            fname = f"{class_name}_export{i}"
            if not args:
                lines.append(
                    f"function {class_name}()\n    ptr = ccall((:{fname}, libfly), Ptr{{Cvoid}}, ())\n    return {class_name}(ptr)\nend\n\n"
                )
            elif args == [("string", None)]:
                lines.append(
                    f"function {class_name}(filename::String)\n    ptr = ccall((:{fname}, libfly), Ptr{{Cvoid}}, (Cstring,), filename)\n    return {class_name}(ptr)\nend\n\n"
                )

        for i, method in enumerate(methods):
            args = method["arg_types"]
            ret = method["ret_type"]

            if args and args[0][0] == "ptr":
                args = args[1:]

            jl_args = []
            ccall_types = ["Ptr{Cvoid}"]
            ccall_vals = ["self.ptr"]
            pre = []

            # Default logic
            min_args = 0
            for j, (_, default) in enumerate(args):
                if default is not None:
                    min_args = j
                    for _, d in args[j:]:
                        if d is None:
                            raise ValueError("Non-default argument after default one")
                    break
            else:
                min_args = len(args)

            for j, (t, default) in enumerate(args):
                var = f"arg{j}"
                if t == "float":
                    jl_args.append(f"{var}" + (f"={default}" if default else ""))
                    pre.append(f"    new{var} = Float32({var})")
                    ccall_types.append("Float32")
                    ccall_vals.append(f"new{var}")
                elif t == "int":
                    jl_args.append(f"{var}::Int" + (f"={default}" if default else ""))
                    ccall_types.append("Cint")
                    ccall_vals.append(var)
                elif t == "string":
                    jl_args.append(
                        f"{var}::String" + (f'="{default}"' if default else "")
                    )
                    ccall_types.append("Cstring")
                    ccall_vals.append(var)
                elif t == "Vec":
                    jl_args.append(f"{var}::Vector{{Float64}}")
                    pre.append(f"    new{var} = Float32.({var})")
                    pre.append(f"    len{var} = length({var})")
                    ccall_types.extend(["Ptr{Float32}", "Cint"])
                    ccall_vals.extend([f"new{var}", f"len{var}"])

            if ret == "complex<float>":
                pre.append("    real = Ref{Float32}()")
                pre.append("    imag = Ref{Float32}()")
                ccall_types += ["Ptr{Float32}", "Ptr{Float32}"]
                ccall_vals += ["real", "imag"]

            # Call operator override
            lines.append(
                f"function (self::{class_name})({', '.join(jl_args)})::{type_map[ret] if ret != 'complex<float>' else 'ComplexF32'}\n"
            )
            lines.extend([line + "\n" for line in pre])
            if ret == "complex<float>":
                lines.append(
                    f"    ccall((:{class_name}_operator_export{i}, libfly), Nothing, ({', '.join(ccall_types)}), {', '.join(ccall_vals)})\n"
                )
                lines.append("    return complex(real[], imag[])\n")
            else:
                lines.append(
                    f"    return ccall((:{class_name}_operator_export{i}, libfly), {type_map[ret]}, ({', '.join(ccall_types)}), {', '.join(ccall_vals)})\n"
                )
            lines.append("end\n\n")

        lines.append(
            f"function Base.finalize(obj::{class_name})\n    destroy!(obj)\nend\n\n"
        )

    # === Main Flow ===
    functions, classes = split_inputs(INPUTS)

    lines.append(f"export {', '.join(sorted(set(exports)))}\n\n")

    for f in functions:
        emit_free_function(f)
    for cname, defs in classes.items():
        emit_class(cname, defs["constructors"], defs["methods"])

    return lines
