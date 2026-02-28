from pathlib import Path
import shutil

def get_extension(language):
    extensions = {
        "python": "py",
        "julia": "jl",
        "c++": "cpp",
        "fortran": "f90",
    }
    return extensions.get(language, "txt")

def make_run_file(language, dir, template_dir):
    extension = get_extension(language)
    file_path = Path(dir) / f"run.{extension}"
    if file_path.exists():
        return
    if language == "python":
        src = template_dir / "python_run.txt"
        shutil.copy2(src, file_path)
    elif language == "julia":
        src = template_dir / "julia_run.txt"
        shutil.copy2(src, file_path)
    elif language == "c++":
        # Determine project root dynamically
        script_dir = Path(__file__).parent  # .../src/config/write
        project_root = script_dir.parent.parent.parent  # .../FFirefly
        config_path = project_root / "build/bin/input.cfg"

        # Read template and replace placeholder
        src = template_dir / "cpp_run.txt"
        with open(src, 'r') as f:
            template = f.read()
        content = template.replace("{CONFIG_PATH}", str(config_path))

        # Write customized file
        with open(file_path, 'w') as f:
            f.write(content)

        # Copy hpp file as-is
        src = template_dir / "hpp_run.txt"
        shutil.copy2(src, Path(dir) / "run.hpp")
    else:
        print("Don't have fortran yet")

def make_test_file(language, dir, template_dir, category=None, calc_name=None, method_name=None):
    extension = get_extension(language)
    file_path = Path(dir) / "tests" / f"test.{extension}"
    if file_path.exists():
        return
    if language == "python":
        src = template_dir / "python_test.txt"
        # Read template and replace placeholders
        with open(src, 'r') as f:
            template = f.read()
        if category and calc_name and method_name:
            content = template.replace("{CATEGORY}", category)
            content = content.replace("{CALCULATION}", calc_name)
            content = content.replace("{METHOD}", method_name)
        else:
            content = template
        with open(file_path, 'w') as f:
            f.write(content)
    elif language == "julia":
        src = template_dir / "julia_test.txt"
        shutil.copy2(src, file_path)
    elif language == "c++":
        src = template_dir / "cpp_test.txt"
        shutil.copy2(src, file_path)
        src = template_dir / "hpp_test.txt"
        shutil.copy2(src, Path(dir) / "tests" / "test.hpp")
    else:
        print("Don't have fortran yet")

def make_readme_file(category, calculation, method, dir, template_dir):
    readme_path = Path(dir) / "README.md"
    if readme_path.exists():
        return

    # Read the template
    template_path = template_dir / "README_template.txt"
    with open(template_path, 'r') as f:
        template = f.read()

    # Replace placeholders
    content = template.replace("{CATEGORY}", category)
    content = content.replace("{CALCULATION}", calculation)
    content = content.replace("{METHOD}", method)

    with open(readme_path, 'w') as f:
        f.write(content)


def make_node_files(category, calculations, category_dir):
    # Build the if-else chain for all calculation/method combinations
    conditions = []
    for calc_name, methods in calculations.items():
        for method_name, language in methods.items():
            if language == "c++":
                conditions.append(
                    f'    if (calculation == "{calc_name}" && method == "{method_name}") {{\n'
                    f'        if (debug) run_cpp_test("{method_name}");\n'
                    f'        else run_cpp_method("{method_name}");\n'
                    f'    }}'
                )
            elif language == "python":
                conditions.append(
                    f'    if (calculation == "{calc_name}" && method == "{method_name}") {{\n'
                    f'        if (debug) run_python_test("{method_name}");\n'
                    f'        else run_python_method("{method_name}");\n'
                    f'    }}'
                )
            elif language == "julia":
                conditions.append(
                    f'    if (calculation == "{calc_name}" && method == "{method_name}") {{\n'
                    f'        if (debug) run_julia_test("{method_name}");\n'
                    f'        else run_julia_method("{method_name}");\n'
                    f'    }}'
                )

    # Join with "else " prefix for all but the first
    if_chain = conditions[0] if conditions else ""
    for cond in conditions[1:]:
        if_chain += "\n    else " + cond.strip()

    node_cpp_path = Path(category_dir) / "node.cpp"
    contents = f"""#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void {category}_wrapper() {{
    printv("Running {category}_wrapper\\n");
{if_chain}
    else {{
        printf("In {category} category, calculation `%s` with method `%s` not recognized\\n", calculation.c_str(), method.c_str());
    }}
}}"""
    with open(node_cpp_path, "w") as f:
        f.write(contents)

    node_hpp_path = Path(category_dir) / "node.hpp"
    contents = f"""#pragma once

    #ifdef __cplusplus
    extern "C" {{
    #endif
    
    void {category}_wrapper();
    
    #ifdef __cplusplus
    }}
    #endif
    """
    with open(node_hpp_path, "w") as f:
        f.write(contents)



def make_folders(categories, base_dir=None):
    """
    Create folder structure for categories.

    Args:
        categories: Dict of {category: {calculation: {method: language}}}
        base_dir: Base directory to create folders in (defaults to ../../)
    """
    if base_dir is None:
        # Get the directory where this script is located
        script_dir = Path(__file__).parent
        base_dir = script_dir / "../../"

    base_dir = Path(base_dir).resolve()
    template_dir = Path(__file__).parent / "templates/" # Where the template .txt files are

    records = []
    for category, calculations in categories.items():
        category_dir = base_dir / category
        category_dir.mkdir(parents=True, exist_ok=True)

        # Create shared directory for this category
        (category_dir / "shared").mkdir(parents=True, exist_ok=True)

        for calc_name, methods in calculations.items():
            calc_dir = category_dir / calc_name
            calc_dir.mkdir(parents=True, exist_ok=True)

            for method_name, language in methods.items():
                method_dir = calc_dir / method_name
                method_dir.mkdir(parents=True, exist_ok=True)
                (method_dir / "tests").mkdir(parents=True, exist_ok=True)

                make_run_file(language, method_dir, template_dir)
                make_test_file(language, method_dir, template_dir, category, calc_name, method_name)
                make_readme_file(category, calc_name, method_name, method_dir, template_dir)

                records.append({
                    'category': category,
                    'calculation': calc_name,
                    'method': method_name,
                    'language': language,
                    'path': str(method_dir)
                })

        # Create node files for this category
        make_node_files(category, calculations, category_dir)

    return records

def update_main_c(categories, base_dir=None):
    """
    Update main.c with category includes and dispatch logic.

    Args:
        categories: Dict of {category: {calculation: {method: language}}}
        base_dir: Base directory (defaults to ../../)
    """
    if base_dir is None:
        script_dir = Path(__file__).parent
        base_dir = script_dir / "../../"

    base_dir = Path(base_dir).resolve()
    main_c_path = base_dir / "main.c"

    # Read the current main.c
    with open(main_c_path, 'r') as f:
        lines = f.readlines()

    # Generate category includes
    category_includes = []
    for category in sorted(categories.keys()):
        category_includes.append(f'#include "{category}/node.hpp"\n')

    # Generate wrapper functions
    wrapper_functions = []
    for category in sorted(categories.keys()):
        wrapper_functions.append(f"void {category}() {{\n")
        wrapper_functions.append(f'    printf("Starting {category} Calculation\\n\\n");\n')
        wrapper_functions.append(f"    {category}_wrapper();\n")
        wrapper_functions.append("}\n\n")

    # Generate dispatch chain
    dispatch_lines = []
    dispatch_lines.append("    /*\n")
    dispatch_lines.append("        * ADDING A CATEGORY OCCURS BELOW\n")
    dispatch_lines.append("        * FOLLOW THE PATTERN\n")
    dispatch_lines.append("    */\n")

    sorted_cats = sorted(categories.keys())
    for i, category in enumerate(sorted_cats):
        if i == 0:
            dispatch_lines.append(f'        if (!strcmp(category, "{category}"))\n')
        else:
            dispatch_lines.append(f'        else if (!strcmp(category, "{category}"))\n')
        dispatch_lines.append(f"            {category}();\n")

    # Add special test case
    dispatch_lines.append(f'        else if (!strcmp(category, "test"))\n')
    dispatch_lines.append(f"            test();\n")
    dispatch_lines.append(f'        else\n')
    dispatch_lines.append(f'            printf("Unknown Category %s\\n", category);\n')

    # Find and replace sections
    new_lines = []
    i = 0
    while i < len(lines):
        # Replace category includes section
        if lines[i].strip() == "// Category nodes below":
            new_lines.append(lines[i])
            i += 1
            # Skip old includes until we hit "// Test nodes below"
            while i < len(lines) and lines[i].strip() != "// Test nodes below":
                i += 1
            # Insert new includes
            for include in category_includes:
                new_lines.append(include)
            new_lines.append("\n")
            continue

        # Replace wrapper functions section
        if lines[i].strip() == "// Global category calls":
            new_lines.append(lines[i])
            new_lines.append("\n")
            i += 1
            # Skip old wrapper functions until we hit print_banner_top
            while i < len(lines) and "void print_banner_top()" not in lines[i]:
                i += 1
            # Insert new wrapper functions
            for func_line in wrapper_functions:
                new_lines.append(func_line)
            continue

        # Replace dispatch section - look for the comment block start
        if lines[i].strip().startswith("/*") and i + 1 < len(lines) and "* ADDING A CATEGORY OCCURS BELOW" in lines[i + 1]:
            # Skip the old dispatch block (comment + all if/else statements + unknown category printf)
            while i < len(lines) and 'printf("Unknown Category %s\\n", category);' not in lines[i]:
                i += 1
            i += 1  # Skip the printf line too
            # Insert new dispatch
            for dispatch_line in dispatch_lines:
                new_lines.append(dispatch_line)
            continue

        new_lines.append(lines[i])
        i += 1

    # Write back to main.c
    with open(main_c_path, 'w') as f:
        f.writelines(new_lines)

    print(f"✓ Updated main.c with {len(categories)} categories")

def write(categories):
    records = make_folders(categories)
    print(f"\nCreated {len(records)} method directories:")
    for record in records:
        print(f"  {record['category']}/{record['calculation']}/{record['method']} ({record['language']})")
    print("\nDone!")

    # Update main.c with category integration
    update_main_c(categories)
