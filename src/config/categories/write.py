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
        src = template_dir / "cpp_run.txt"
        shutil.copy2(src, file_path)
        src = template_dir / "hpp_run.txt"
        shutil.copy2(src, Path(dir) / "run.hpp")
    else:
        print("Don't have fortran yet")

def make_test_file(language, dir, template_dir):
    extension = get_extension(language)
    file_path = Path(dir) / "tests" / f"test.{extension}"
    if file_path.exists():
        return
    if language == "python":
        src = template_dir / "python_test.txt"
        shutil.copy2(src, file_path)
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
                conditions.append(f'    if (calculation == "{calc_name}" && method == "{method_name}") run_cpp_method("{method_name}");')
            elif language == "python":
                conditions.append(f'    if (calculation == "{calc_name}" && method == "{method_name}") run_python_method("{method_name}");')
            elif language == "julia":
                conditions.append(f'    if (calculation == "{calc_name}" && method == "{method_name}") run_julia_method("{method_name}");')

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
                make_test_file(language, method_dir, template_dir)
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

def write(categories):
    records = make_folders(categories)
    print(f"\nCreated {len(records)} method directories:")
    for record in records:
        print(f"  {record['category']}/{record['calculation']}/{record['method']} ({record['language']})")
    print("\nDone!")
