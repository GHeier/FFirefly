import re
from collections import defaultdict
from write import write_export, write_import


def parse_function(line):
    # Match return type, function name, and arguments
    func_pattern = (
        r"([a-zA-Z_:<>0-9]+(?:\s*\*?)?)\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(([^)]*)\)"
    )
    match = re.match(func_pattern, line.strip())
    if match:
        return_type = match.group(1).strip()
        func_name = match.group(2).strip()
        args = match.group(3).strip()
        arg_types = [
            arg.strip().rsplit(" ", 1)[0] for arg in args.split(",") if arg.strip()
        ]
        return func_name, {"args": arg_types, "return": return_type}
    return None


def parse_class_variables(class_body):
    # Only get lines that look like variable declarations
    variables = {}
    lines = class_body.splitlines()
    for line in lines:
        line = line.strip().strip(";")
        if not line or "(" in line or ")" in line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            var_type = " ".join(parts[:-1])
            var_name = parts[-1]
            variables[var_name] = var_type
    return variables


def parse_hpp_file(filename):
    functions = {}
    classes = {}

    with open(filename, "r") as f:
        content = f.read()

    # Remove comments
    content = re.sub(r"//.*?$|/\*.*?\*/", "", content, flags=re.DOTALL | re.MULTILINE)

    # Extract class definitions
    class_pattern = re.compile(r"class\s+(\w+)\s*{([^}]*)};", re.DOTALL)
    for class_match in class_pattern.finditer(content):
        class_name = class_match.group(1)
        class_body = class_match.group(2)
        classes[class_name] = parse_class_variables(class_body)

    # Extract free-standing functions
    lines = content.splitlines()
    for line in lines:
        if "(" in line and ")" in line and line.strip().endswith(";"):
            parsed = parse_function(line)
            if parsed:
                name, info = parsed
                functions[name] = info

    return functions, classes


FILES = {"hamiltonian/band_structure.hpp"}

for file in FILES:
    filename = "../" + file
    functions, classes = parse_hpp_file(filename)
    write_export.write_export_funcs(functions, FILES)
# print("Functions:")
# for name, info in functions.items():
#    print(f"{name}: {info}")

# print("\nClasses:")
# for name, vars in classes.items():
#    print(f"{name}: {vars}")
