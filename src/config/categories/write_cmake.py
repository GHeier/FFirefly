"""
CMake code generator for FFirefly categories.

Generates CMakeLists.txt sections for:
1. base.so shared library source lists
2. C++ method executable targets
"""

from pathlib import Path
import sys

# Import CATEGORIES from parent module if run as part of categories.py
# Otherwise it will be passed as a parameter
try:
    import __main__
    if hasattr(__main__, 'CATEGORIES'):
        CATEGORIES = __main__.CATEGORIES
    else:
        CATEGORIES = None
except:
    CATEGORIES = None


def get_base_library_sources():
    """
    Generate list of source files for base.so library.

    Includes:
    - algorithms/ (all files)
    - config/ (all files)
    - hamiltonian/models/ (only model files, not calculation implementations)
    - module/exports/ (C++ export wrappers)
    - objects/ (all files)
    """
    base_dirs = [
        "src/algorithms",
        "src/config",
        "src/hamiltonian",  # Include all of hamiltonian except method implementations
        "src/module/exports",
        "src/objects",
    ]

    cmake_code = []
    cmake_code.append("# Base library sources")
    cmake_code.append("set(BASE_SOURCES")

    # Get project root (3 levels up from this file: categories/ -> config/ -> src/ -> root)
    project_root = Path(__file__).parent.parent.parent.parent
    project_root = project_root.resolve()

    found_files = 0
    for base_path in base_dirs:
        dir_path = project_root / base_path
        if not dir_path.exists():
            print(f"Warning: Base directory not found: {dir_path}")
            continue

        # Find all C/C++ source files, excluding archive directories and Fortran files
        for ext in ["*.c", "*.cpp"]:
            for file in dir_path.rglob(ext):
                # Skip Fortran files (they're in a separate library)
                if file.suffix == '.f90':
                    continue
                # Skip archive directories
                if "/archive/" in str(file) or "\\archive\\" in str(file):
                    continue

                # Skip method implementation files (they're separate executables)
                # Match pattern: src/hamiltonian/{calculation}/{method}/
                file_str = str(file)
                is_method_file = False
                if base_path == "src/hamiltonian":
                    # Check if this is in a method directory (3+ levels deep)
                    rel_to_base = file.relative_to(dir_path)
                    parts = rel_to_base.parts
                    # If path is like DOS/tetrahedra/run.cpp or FS/tetrahedra/fs.cpp
                    if len(parts) >= 3 and parts[0] in ["DOS", "FS", "generate"]:
                        is_method_file = True

                if is_method_file:
                    continue

                # Get relative path from project root
                rel_path = file.relative_to(project_root)
                cmake_code.append(f"    {rel_path}")
                found_files += 1

    cmake_code.append(")")
    cmake_code.append("")

    print(f"   Found {found_files} source files for base.so")
    return "\n".join(cmake_code)


def get_cpp_method_executables(categories):
    """
    Generate CMake executable targets for all C++ methods.

    For each C++ method in categories, creates:
    - An executable target: {category}_{calculation}_{method}.exe
    - Links against base.so
    - Includes all .c/.cpp files in the method directory
    """
    cmake_code = []
    cmake_code.append("# C++ method executables")
    cmake_code.append("")

    # Get project root
    project_root = Path(__file__).parent.parent.parent.parent
    project_root = project_root.resolve()
    src_dir = project_root / "src"

    cpp_method_count = 0
    for category, calculations in categories.items():
        for calculation, methods in calculations.items():
            for method, language in methods.items():
                if language != "c++":
                    continue

                # Executable name
                exe_name = f"{category}_{calculation}_{method}"

                # Method directory
                method_dir = src_dir / category / calculation / method
                if not method_dir.exists():
                    print(f"   Warning: Method directory not found: {method_dir}")
                    continue

                # Find all source files in method directory
                sources = []
                for ext in ["*.c", "*.cpp"]:
                    for file in method_dir.rglob(ext):
                        # Skip test files
                        if "/tests/" in str(file) or "\\tests\\" in str(file):
                            continue
                        rel_path = file.relative_to(project_root)
                        sources.append(str(rel_path))

                if not sources:
                    print(f"   Warning: No source files found for {exe_name}")
                    continue

                # Generate CMake target
                cmake_code.append(f"# {category}/{calculation}/{method}")
                cmake_code.append(f"add_executable({exe_name}")
                for src in sources:
                    cmake_code.append(f"    {src}")
                cmake_code.append(")")
                cmake_code.append(f"target_link_libraries({exe_name} PRIVATE base)")
                cmake_code.append(f"set_target_properties({exe_name} PROPERTIES")
                cmake_code.append(f"    OUTPUT_NAME \"{exe_name}.exe\"")
                cmake_code.append(f"    RUNTIME_OUTPUT_DIRECTORY \"${{CMAKE_BINARY_DIR}}/bin\"")
                cmake_code.append(f"    BUILD_RPATH \"${{CMAKE_LIBRARY_OUTPUT_DIRECTORY}}\"")
                cmake_code.append(f"    INSTALL_RPATH \"${{CMAKE_LIBRARY_OUTPUT_DIRECTORY}}\"")
                cmake_code.append(")")
                cmake_code.append("")
                cpp_method_count += 1

    print(f"   Generated {cpp_method_count} C++ method executables")
    return "\n".join(cmake_code)


def update_cmakelists(categories):
    """
    Update CMakeLists.txt with auto-generated sections.

    Replaces content between marker comments:
    - # BEGIN AUTO-GENERATED BASE LIBRARY SOURCES
    - # END AUTO-GENERATED BASE LIBRARY SOURCES
    - # BEGIN AUTO-GENERATED CPP METHOD EXECUTABLES
    - # END AUTO-GENERATED CPP METHOD EXECUTABLES
    """
    cmake_path = Path(__file__).parent / "../../../CMakeLists.txt"
    cmake_path = cmake_path.resolve()

    if not cmake_path.exists():
        print(f"Error: CMakeLists.txt not found at {cmake_path}")
        return False

    with open(cmake_path, 'r') as f:
        content = f.read()

    # Generate new sections
    base_sources = get_base_library_sources()
    cpp_executables = get_cpp_method_executables(categories)

    # Replace base library sources section
    base_marker_start = "# BEGIN AUTO-GENERATED BASE LIBRARY SOURCES"
    base_marker_end = "# END AUTO-GENERATED BASE LIBRARY SOURCES"

    if base_marker_start in content and base_marker_end in content:
        start_idx = content.find(base_marker_start)
        end_idx = content.find(base_marker_end) + len(base_marker_end)

        new_section = f"{base_marker_start}\n{base_sources}\n{base_marker_end}"
        content = content[:start_idx] + new_section + content[end_idx:]
    else:
        print(f"Warning: Base library source markers not found in CMakeLists.txt")

    # Replace C++ executables section
    exe_marker_start = "# BEGIN AUTO-GENERATED CPP METHOD EXECUTABLES"
    exe_marker_end = "# END AUTO-GENERATED CPP METHOD EXECUTABLES"

    if exe_marker_start in content and exe_marker_end in content:
        start_idx = content.find(exe_marker_start)
        end_idx = content.find(exe_marker_end) + len(exe_marker_end)

        new_section = f"{exe_marker_start}\n{cpp_executables}\n{exe_marker_end}"
        content = content[:start_idx] + new_section + content[end_idx:]
    else:
        print(f"Warning: C++ executable markers not found in CMakeLists.txt")

    # Write back
    with open(cmake_path, 'w') as f:
        f.write(content)

    print(f"✓ Updated CMakeLists.txt")
    return True


if __name__ == "__main__":
    print("Generating CMake configuration from categories...")
    update_cmakelists(CATEGORIES)
    print("Done!")
