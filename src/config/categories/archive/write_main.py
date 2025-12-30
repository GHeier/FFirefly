#!/usr/bin/env python3
"""
Generates main.c sections for category integration.
This script creates the code snippets needed to integrate categories into main.c
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
CATEGORIES_WITH_TESTS = categories_module.CATEGORIES_WITH_TESTS


def generate_includes():
    """Generate #include statements for all category nodes"""
    includes = []
    includes.append("// Category nodes below")

    for category in sorted(CATEGORIES.keys()):
        includes.append(f'#include "{category}/node.hpp"')

    return "\n".join(includes)


def generate_test_includes():
    """Generate #include statements for all test nodes"""
    includes = []
    includes.append("\n// Test nodes below")

    for category in CATEGORIES_WITH_TESTS:
        includes.append(f'#include "{category}/tests/all.hpp"')

    return "\n".join(includes)


def generate_function_definitions():
    """Generate function definitions that call category wrappers"""
    functions = []
    functions.append("\n// Global category calls\n")

    for category in sorted(CATEGORIES.keys()):
        func_name = category.replace("/", "_")
        functions.append(f"""void {func_name}() {{
    printf("Starting {category.title()} Calculation\\n\\n");
    {func_name}_wrapper();
}}
""")

    return "\n".join(functions)


def generate_dispatcher():
    """Generate the if-else dispatcher for categories"""
    lines = []
    lines.append("    /*")
    lines.append("        * ADDING A CATEGORY OCCURS BELOW")
    lines.append("        * FOLLOW THE PATTERN")
    lines.append("    */")

    for i, category in enumerate(sorted(CATEGORIES.keys())):
        func_name = category.replace("/", "_")
        if i == 0:
            lines.append(f'        if (!strcmp(category, "{category}"))')
        else:
            lines.append(f'        else if (!strcmp(category, "{category}"))')
        lines.append(f'            {func_name}();')

    lines.append('        else if (!strcmp(category, "test"))')
    lines.append('            test();')
    lines.append('        else')
    lines.append('            printf("Unknown Category %s\\n", category);')

    return "\n".join(lines)


def generate_test_function():
    """Generate the test() function with all test calls"""
    lines = []
    lines.append("// Global test calls")
    lines.append("void test() {")
    lines.append('    printf("Starting Test Calculations\\n");')
    lines.append(f"\n    int num_tests = {len(CATEGORIES_WITH_TESTS)};")
    lines.append("\n    bool all_tests[num_tests];")

    for i, category in enumerate(CATEGORIES_WITH_TESTS):
        # Handle special naming cases
        if category == "objects":
            test_func = "object_tests"
        elif category == "algorithms":
            test_func = "algorithm_tests"
        else:
            test_func = category.replace("/", "_") + "_tests"
        lines.append(f"    all_tests[{i}] = {test_func}();")

    lines.append('\n    printf("\\n");')
    lines.append('    print_test_results(all_tests, num_tests, "Test Categories");')
    lines.append("}")

    return "\n".join(lines)


def generate_main_sections():
    """Generate all main.c sections and write to output file"""
    output_dir = os.path.join(os.path.dirname(__file__), "archive")
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "main_sections.txt")

    with open(output_file, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("GENERATED CODE FOR main.c\n")
        f.write("=" * 70 + "\n\n")

        f.write("// SECTION 1: Includes (add after config includes)\n")
        f.write(generate_includes() + "\n")
        f.write(generate_test_includes() + "\n\n")

        f.write("=" * 70 + "\n\n")
        f.write("// SECTION 2: Function definitions (before main)\n")
        f.write(generate_function_definitions() + "\n")
        f.write(generate_test_function() + "\n\n")

        f.write("=" * 70 + "\n\n")
        f.write("// SECTION 3: Dispatcher (inside main loop)\n")
        f.write(generate_dispatcher() + "\n\n")

        f.write("=" * 70 + "\n")

    print(f"Main sections generated: {output_file}")
    return output_file


if __name__ == "__main__":
    generate_main_sections()
