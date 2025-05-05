def create_julia_function(INPUTS):
    lines = []
    exports = []

    def parse_arg(arg):
        if "=" in arg:
            t, default = arg.split("=")
            return t.strip(), default.strip()
        return arg.strip(), None

    type_map = {
        "int": "Cint",
        "float": "Cfloat",
        "double": "Cdouble",
        "string": "Cstring",
        "ptr": "Ptr{Cvoid}",
        "complex<float>": "Nothing",  # returned via Ref
        "Vec": ("Ptr{Cfloat}", "Cint"),
        "None": "Nothing",
    }

    julia_types = {
        "int": "Int",
        "float": "Float32",
        "double": "Float64",
        "string": "String",
        "ptr": "Ptr{Cvoid}",
        "complex<float>": "ComplexF32",
        "Vec": "Vector{Float64}",
    }

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

    lines.append(f"export {', '.join(sorted(exports))}\n\n")

    # Free functions
    for entry in functions:
        name = entry["name"]
        ret_type = entry["ret_type"]
        arg_types = entry["arg_types"]

        jl_args = []
        ccall_args = []
        ccall_types = []
        pre = []

        for i, (t, default) in enumerate(arg_types):
            var = f"arg{i}"
            if t == "Vec":
                jl_args.append(f"{var}::Vector{{Float64}}")
                pre.append(f"    new{var} = Float32.({var})")
                ccall_types.extend(["Ptr{Cfloat}", "Cint"])
                ccall_args.extend([f"new{var}", f"length({var})"])
            else:
                base = julia_types[t]
                jl = f"{var}::{base}"
                if default:
                    jl += f" = {default}"
                jl_args.append(jl)
                ccall_types.append(type_map[t])
                ccall_args.append(var)

        if ret_type == "complex<float>":
            lines.append(f"function {name}({', '.join(jl_args)})::ComplexF32\n")
            lines.extend([p + "\n" for p in pre])
            lines.append("    re = Ref{Cfloat}()\n")
            lines.append("    im = Ref{Cfloat}()\n")
            lines.append(
                f"    ccall((:{name}_export0, libfly), Nothing, (Ref{{Cfloat}}, Ref{{Cfloat}}, {', '.join(ccall_types)}), re, im, {', '.join(ccall_args)})\n"
            )
            lines.append("    return ComplexF32(re[], im[])\n")
            lines.append("end\n\n")
        else:
            jl_ret = type_map[ret_type]
            lines.append(f"function {name}({', '.join(jl_args)})::{jl_ret}\n")
            lines.extend([p + "\n" for p in pre])
            lines.append(
                f"    return ccall((:{name}_export0, libfly), {jl_ret}, ({', '.join(ccall_types)}), {', '.join(ccall_args)})\n"
            )
            lines.append("end\n\n")

    # Classes
    for cls, defs in classes.items():
        lines.append(f"mutable struct {cls}\n")
        lines.append("    ptr::Ptr{Cvoid}\n")
        lines.append("end\n\n")

        # Constructors
        for i, ctor in enumerate(defs["constructors"]):
            args = ctor["arg_types"]
            cname = f"{cls}_export{i}"
            if not args:
                lines.append(f"function {cls}()\n")
                lines.append(f"    ptr = ccall((:{cname}, libfly), Ptr{{Cvoid}}, ())\n")
                lines.append(f"    return {cls}(ptr)\n")
                lines.append("end\n\n")
            elif args == [("string", None)]:
                lines.append(f"function {cls}(filename::String)\n")
                lines.append(
                    f"    ptr = ccall((:{cname}, libfly), Ptr{{Cvoid}}, (Cstring,), filename)\n"
                )
                lines.append(f"    return {cls}(ptr)\n")
                lines.append("end\n\n")

        # Operator overloads (call-style)
        for i, method in enumerate(defs["methods"]):
            ret_type = method["ret_type"]
            method_args = method["arg_types"]
            if method_args and method_args[0][0] == "ptr":
                method_args = method_args[1:]

            jl_args = [f"self::{cls}"]
            ccall_args = ["self.ptr"]
            ccall_types = ["Ptr{Cvoid}"]
            pre = []

            for j, (t, default) in enumerate(method_args):
                var = f"arg{j}"
                if t == "float":
                    jl = f"{var}::Float32"
                    if default:
                        jl += f" = {default}"
                    jl_args.append(jl)
                    pre.append(f"    new{var} = Cfloat({var})")
                    ccall_types.append("Cfloat")
                    ccall_args.append(f"new{var}")
                elif t == "int":
                    jl_args.append(f"{var}::Int")
                    ccall_types.append("Cint")
                    ccall_args.append(var)
                elif t == "Vec":
                    jl_args.append(f"{var}::Vector{{Float64}}")
                    pre.append(f"    new{var} = Float32.({var})")
                    pre.append(f"    len = length({var})")
                    ccall_types.extend(["Ptr{Cfloat}", "Cint"])
                    ccall_args.extend([f"new{var}", "len"])

            mname = f"{cls}_operator_export{i}"

            if ret_type == "complex<float>":
                lines.append(
                    f"function ({jl_args[0]})({', '.join(jl_args[1:])})::ComplexF32\n"
                )
                lines.extend([p + "\n" for p in pre])
                lines.append("    re = Ref{Cfloat}()\n")
                lines.append("    im = Ref{Cfloat}()\n")
                lines.append(
                    f"    ccall((:{mname}, libfly), Nothing, (Ref{{Cfloat}}, Ref{{Cfloat}}, {', '.join(ccall_types)}), re, im, {', '.join(ccall_args)})\n"
                )
                lines.append("    return ComplexF32(re[], im[])\n")
                lines.append("end\n\n")
            else:
                jl_ret = type_map[ret_type]
                lines.append(
                    f"function ({jl_args[0]})({', '.join(jl_args[1:])})::{jl_ret}\n"
                )
                lines.extend([p + "\n" for p in pre])
                lines.append(
                    f"    return ccall((:{mname}, libfly), {jl_ret}, ({', '.join(ccall_types)}), {', '.join(ccall_args)})\n"
                )
                lines.append("end\n\n")

        lines.append(f"function Base.finalize(obj::{cls})\n")
        lines.append("    destroy!(obj)\n")
        lines.append("end\n\n")

    return lines
