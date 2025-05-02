def create_julia_function(INPUTS):
    lines = []
    exports = []

    type_map = {
        "int": "Cint",
        "float": "Float32",
        "double": "Float64",
        "string": "Cstring",
        "ptr": "Ptr{Cvoid}",
        "complex<float>": "NTuple{2, Float32}",
        "None": "Nothing",
        "Vec": ("Ptr{Float32}", "Cint"),
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
                exports.append(base_name)
            elif "_operator" in base_name:
                class_name = base_name.split("_operator")[0]
                if class_name in classes:
                    classes[class_name]["methods"].append(entry)
            else:
                functions.append(entry)
                exports.append(base_name)

    # Export line
    lines.append(f"export {', '.join(sorted(exports))}\n\n")

    # Free functions
    for entry in functions:
        name = entry["name"]
        ret_type = type_map[entry["ret_type"]]
        arg_types = entry["arg_types"]

        julia_args = []
        ccall_types = []
        ccall_vals = []
        pre = []

        for i, t in enumerate(arg_types):
            var = f"arg{i}"
            if t == "Vec":
                julia_args.append(f"{var}::Vector{{Float64}}")
                pre.append(f"    new{var} = Float32.({var})")
                ccall_types.extend(["Ptr{Float32}", "Cint"])
                ccall_vals.extend([f"new{var}", f"length({var})"])
            else:
                julia_args.append(f"{var}::{julia_types[t]}")
                ccall_types.append(type_map[t])
                ccall_vals.append(var)

        lines.append(f"function {name}({', '.join(julia_args)})\n")
        lines.extend([p + "\n" for p in pre])
        ret_str = f"{ret_type}" if ret_type != "Nothing" else "Nothing"
        call = f'return ccall((:{name}_export0, libfly), {ret_str}, ({", ".join(ccall_types)}), {", ".join(ccall_vals)})'
        lines.append(f"    {call}\n")
        lines.append("end\n\n")

    # Object wrappers
    for cls, defn in classes.items():
        lines.append(f"mutable struct {cls}\n")
        lines.append("    ptr::Ptr{Cvoid}\n")
        lines.append("end\n\n")

        for i, ctor in enumerate(defn["constructors"]):
            args = ctor["arg_types"]
            cname = f"{cls}_export{i}"
            if not args:
                lines.append(f"function {cls}()\n")
                lines.append(f"    ptr = ccall((:{cname}, libfly), Ptr{{Cvoid}}, ())\n")
                lines.append(f"    return {cls}(ptr)\n")
                lines.append("end\n\n")
            elif args == ["string"]:
                lines.append(f"function {cls}(filename::String)\n")
                lines.append(
                    f"    ptr = ccall((:{cname}, libfly), Ptr{{Cvoid}}, (Cstring,), filename)\n"
                )
                lines.append(f"    return {cls}(ptr)\n")
                lines.append("end\n\n")

        for i, method in enumerate(defn["methods"]):
            args = method["arg_types"]
            mname = f"{cls}_operator_export{i}"

            julia_sig = [f"self::{cls}"]
            c_types = ["Ptr{Cvoid}"]
            c_vals = ["self.ptr"]
            pre = []

            for j, t in enumerate(args):
                var = f"arg{j}"
                if t == "float":
                    pre.append(f"    new{var} = Float32({var})")
                    c_types.append("Float32")
                    c_vals.append(f"new{var}")
                    julia_sig.append(f"{var}")
                elif t == "int":
                    c_types.append("Cint")
                    c_vals.append(var)
                    julia_sig.append(f"{var}")
                elif t == "Vec":
                    pre.append(f"    new{var} = Float32.({var})")
                    pre.append(f"    len = length({var})")
                    c_types.extend(["Ptr{Float32}", "Cint"])
                    c_vals.extend([f"new{var}", "len"])
                    julia_sig.append(f"{var}::Vector{{Float64}}")
                else:
                    julia_sig.append(f"{var}")

            ret_jl = type_map[method["ret_type"]]
            lines.append(f"function ({', '.join(julia_sig)})::{ret_jl}\n")
            lines.extend([p + "\n" for p in pre])
            lines.append(
                f"    return ccall((:{mname}, libfly), {ret_jl}, ({', '.join(c_types)}), {', '.join(c_vals)})\n"
            )
            lines.append("end\n\n")

        lines.append(f"function Base.finalize(obj::{cls})\n")
        lines.append(f"    destroy!(obj)\n")
        lines.append("end\n\n")

    return lines
