#!/usr/bin/env python3
"""
Generates node.cpp and node.hpp files for all categories.
These files serve as wrappers that dispatch to the appropriate calculation functions.
"""

import os
import sys

# Load categories from parent module
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("categories_def",
                                                os.path.join(os.path.dirname(os.path.dirname(__file__)), "categories.py"))
categories_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(categories_module)
CATEGORIES = categories_module.CATEGORIES


def get_file_extension(filename):
    """Extract file extension from filename"""
    return os.path.splitext(filename)[1]


def get_function_name(filename):
    """Extract function name from filename (without extension)"""
    return os.path.splitext(filename)[0]


def generate_node_hpp(category):
    """Generate node.hpp header file for a category"""
    func_name = category.replace("/", "_")

    content = f"""#pragma once

#ifdef __cplusplus
extern "C" {{
#endif

void {func_name}_wrapper();

#ifdef __cplusplus
}}
#endif
"""
    return content


def generate_node_cpp(category, calculations):
    """Generate node.cpp implementation file for a category"""
    func_name = category.replace("/", "_")

    # Determine what to include based on file types used
    includes = []
    includes.append('#include "node.hpp"')
    includes.append('#include "../config/load/cpp_config.hpp"')

    has_python = False
    has_julia = False
    cpp_files = set()

    for calc_name, methods in calculations.items():
        for method_name, filename in methods.items():
            ext = get_file_extension(filename)
            if ext == ".py":
                has_python = True
            elif ext == ".jl":
                has_julia = True
            elif ext == ".cpp":
                # Use the actual filename for the header, not calc_name
                header = get_function_name(filename) + ".hpp"
                cpp_files.add(header)

    if has_python:
        includes.append('#include "../config/load/py_interface.h"')
    if has_julia:
        includes.append('#include "../config/load/jl_interface.h"')

    for header in sorted(cpp_files):
        includes.append(f'#include "{header}"')

    includes.append("#include <iostream>")
    includes.append("\nusing namespace std;")

    # Generate wrapper functions for Python/Julia calculations
    helper_functions = []

    for calc_name, methods in sorted(calculations.items()):
        for method_name, filename in sorted(methods.items()):
            ext = get_file_extension(filename)
            if ext in [".py", ".jl"]:
                func = get_function_name(filename)
                if method_name is None:
                    helper_func_name = f"{calc_name}_impl"
                else:
                    helper_func_name = f"{calc_name}_{method_name}_impl"

                if ext == ".jl":
                    # Julia function
                    module_name = "".join(word.title() for word in func.split("_"))
                    helper_functions.append(f"""
void {helper_func_name}() {{
    string folder = "{category}/";
    string filename = "{func}";
    string module = "{module_name}";
    string function = "main";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}}""")
                else:
                    # Python function
                    helper_functions.append(f"""
void {helper_func_name}() {{
    string folder = "{category}";
    string filename = "{func}";
    string function = "main";
    call_python_func(folder.c_str(), filename.c_str(), function.c_str());
}}""")

    # Generate main wrapper function with dispatcher
    dispatcher_lines = []

    for i, (calc_name, methods) in enumerate(sorted(calculations.items())):
        prefix = "if" if i == 0 else "else if"

        # Check if this calculation has methods or direct call
        if None in methods:
            # Direct call without method checking
            filename = methods[None]
            ext = get_file_extension(filename)
            # For C++, use calculation name as function name (not filename)
            # This allows multiple calculations to use the same file
            func = calc_name if ext == ".cpp" else get_function_name(filename)

            dispatcher_lines.append(f'    {prefix} (calculation == "{calc_name}") {{')
            if ext == ".cpp":
                dispatcher_lines.append(f'        {func}();')
            else:
                dispatcher_lines.append(f'        {calc_name}_impl();')
            dispatcher_lines.append('    }')
        else:
            # Has methods - need nested if-else for method dispatch
            dispatcher_lines.append(f'    {prefix} (calculation == "{calc_name}") {{')

            for j, (method_name, filename) in enumerate(sorted(methods.items())):
                ext = get_file_extension(filename)
                # For C++, try to use calc_name_method as function name
                if ext == ".cpp":
                    func = f"{calc_name}_{method_name}" if method_name else calc_name
                else:
                    func = get_function_name(filename)
                method_prefix = "if" if j == 0 else "else if"

                dispatcher_lines.append(f'        {method_prefix} (method == "{method_name}") {{')
                if ext == ".cpp":
                    dispatcher_lines.append(f'            {func}();')
                else:
                    dispatcher_lines.append(f'            {calc_name}_{method_name}_impl();')
                dispatcher_lines.append('        }')

            dispatcher_lines.append('        else')
            dispatcher_lines.append(f'            cout << "method \\"" << method << "\\" not recognized for calculation {calc_name}" << endl;')
            dispatcher_lines.append('    }')

    dispatcher_lines.append('    else')
    dispatcher_lines.append('        cout << "calculation \\"" << calculation << "\\" not recognized for category ' + category + '" << endl;')

    wrapper_function = f"""
/**
 * Wrapper function for {category} category
 * Dispatches to appropriate calculation/method based on config variables
 */
extern "C" void {func_name}_wrapper() {{
    printv("Running {func_name}_wrapper\\n");
{chr(10).join(dispatcher_lines)}
}}"""

    # Assemble full content
    content = "\n".join(includes)
    content += "\n" + "\n".join(helper_functions)
    content += "\n" + wrapper_function + "\n"

    return content


def generate_stub_file(category, calc_name, method_name, filename):
    """Generate stub implementation file for a calculation"""
    ext = get_file_extension(filename)
    func_name = get_function_name(filename)

    if ext == ".cpp":
        # C++ stub
        content = f"""#include "{func_name}.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>

using namespace std;

/**
 * {category}/{calc_name}"""
        if method_name:
            content += f" - method: {method_name}"
        content += f"""
 *
 * TODO: Implement this function
 */
void {func_name}() {{
    cout << "Running {func_name}() - NOT YET IMPLEMENTED" << endl;
    // TODO: Add implementation here
}}
"""
        # Also generate header
        header_content = f"""#pragma once

void {func_name}();
"""
        return content, header_content

    elif ext == ".py":
        # Python stub
        content = f"""#!/usr/bin/env python3
\"\"\"
{category}/{calc_name}"""
        if method_name:
            content += f" - method: {method_name}"
        content += f"""

TODO: Implement this module
\"\"\"


def main():
    \"\"\"Main entry point for {func_name}\"\"\"
    print("Running {func_name} - NOT YET IMPLEMENTED")
    # TODO: Add implementation here
    pass


if __name__ == "__main__":
    main()
"""
        return content, None

    elif ext == ".jl":
        # Julia stub
        module_name = "".join(word.title() for word in func_name.split("_"))
        content = f"""# {category}/{calc_name}"""
        if method_name:
            content += f" - method: {method_name}"
        content += f"""
# TODO: Implement this module

module {module_name}

export main

\"\"\"
    main()

Main entry point for {module_name}
TODO: Implement this function
\"\"\"
function main()
    println("Running {module_name}.main() - NOT YET IMPLEMENTED")
    # TODO: Add implementation here
end

end # module
"""
        return content, None

    return None, None


def create_implementation_files():
    """Create stub implementation files for all calculations"""
    src_dir = os.path.join(os.path.dirname(__file__), "..", "..")

    for category, calculations in CATEGORIES.items():
        category_dir = os.path.join(src_dir, category)
        os.makedirs(category_dir, exist_ok=True)

        for calc_name, methods in calculations.items():
            for method_name, filename in methods.items():
                # Check if file already exists
                impl_file = os.path.join(category_dir, filename)

                if os.path.exists(impl_file):
                    print(f"Skipping (exists): {impl_file}")
                    continue

                # Generate stub content
                content, header_content = generate_stub_file(category, calc_name, method_name, filename)

                if content:
                    # Write implementation file
                    with open(impl_file, "w") as f:
                        f.write(content)
                    print(f"Created: {impl_file}")

                    # Write header file if C++
                    if header_content:
                        header_file = os.path.join(category_dir, get_function_name(filename) + ".hpp")
                        if not os.path.exists(header_file):
                            with open(header_file, "w") as f:
                                f.write(header_content)
                            print(f"Created: {header_file}")


def generate_all_nodes():
    """Generate node files for all categories"""
    src_dir = os.path.join(os.path.dirname(__file__), "..", "..")

    for category, calculations in CATEGORIES.items():
        category_dir = os.path.join(src_dir, category)

        # Create output directory for archive files (excluded from build)
        output_dir = os.path.join(os.path.dirname(__file__), "archive", category)
        os.makedirs(output_dir, exist_ok=True)

        # Generate and write node.hpp
        hpp_content = generate_node_hpp(category)
        hpp_file = os.path.join(output_dir, "node.hpp")
        with open(hpp_file, "w") as f:
            f.write(hpp_content)
        print(f"Generated: {hpp_file}")

        # Generate and write node.cpp
        cpp_content = generate_node_cpp(category, calculations)
        cpp_file = os.path.join(output_dir, "node.cpp")
        with open(cpp_file, "w") as f:
            f.write(cpp_content)
        print(f"Generated: {cpp_file}")


if __name__ == "__main__":
    generate_all_nodes()
    create_implementation_files()
