#!/usr/bin/env python3
"""
Category definitions for FFirefly
This file defines all categories, their calculations, methods, and implementation files.

Structure:
CATEGORIES = {
    "category_name": {
        "calculation_name": {
            "method_name": "filename.ext",  # .py, .jl, or .cpp
            # OR if no methods needed:
            None: "filename.ext",  # Direct calculation without method dispatch
        }
    }
}

File extensions determine how calculations are called:
- .py  -> call_python_func()
- .jl  -> call_julia_func()
- .cpp -> direct function call
"""

CATEGORIES = {
    "hamiltonian": {
        "fs": {
            None: "fs.cpp",  # No method needed, direct call
        },
        "dos": {
            "libtetrabz": "dos_tetrabz.cpp",
            "python": "dos_python.py",
        },
        "bands": {
            None: "band_structure.cpp",
        },
    },
    "many_body": {
        "vertex": {
            "FLEX": "vertex.cpp",
        },
        "self_energy": {
            "sparse_ir": "self_energy.cpp",
        },
        "renormalization": {
            None: "renormalization.cpp",  # Has internal method dispatch
        },
        "response": {
            "libtetrabz": "response_tetrabz.cpp",
            "sparse_ir": "response_ir.jl",
        },
        "loop": {
            None: "many_body_loop.jl",
        },
        "triqs": {
            None: "many_body_triqs.py",
        },
    },
    "superconductor": {
        "bcs": {
            None: "superconductor.cpp",  # bcs() function
        },
        "eliashberg": {
            None: "superconductor.cpp",  # eliashberg() function
        },
        "linearized_eliashberg": {
            None: "superconductor.cpp",  # linearized_eliashberg() function
        },
        "debug": {
            None: "superconductor.cpp",  # debug() function
        },
    },
}

# List of categories that have test suites
CATEGORIES_WITH_TESTS = [
    "hamiltonian",
    "objects",
    "algorithms",
    "config/load",
    "module",
]

if __name__ == "__main__":
    # Run code generators when this file is executed
    import os
    import sys
    import importlib.util

    # Load generator modules directly from file paths
    generators_dir = os.path.join(os.path.dirname(__file__), "categories")

    # Load write_main
    spec = importlib.util.spec_from_file_location("write_main",
                                                    os.path.join(generators_dir, "write_main.py"))
    write_main = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(write_main)

    # Load write_nodes
    spec = importlib.util.spec_from_file_location("write_nodes",
                                                    os.path.join(generators_dir, "write_nodes.py"))
    write_nodes = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(write_nodes)

    print("Generating category code...")
    write_main.generate_main_sections()
    write_nodes.generate_all_nodes()
    print("\nCreating stub implementation files...")
    write_nodes.create_implementation_files()
    print("\nDone! Category code generated successfully.")
