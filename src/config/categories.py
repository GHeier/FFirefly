# Setup as:
# CATEGORIES = {
#   category1 = {
#       calculation1 = {
#           method1: "python",
#           method2: "julia",
#           method3: "c++"
#       },
#       calculation2 = ...
#   },
#   category2 = ...

CATEGORIES = {
    "hamiltonian": { # CATEGORY NAME
        "DOS": { # CALCULATION NAME
            "gaussian": "python", # METHOD NAME: IMPLEMENTATION LANGUAGE
            "tetrahedra": "c++"
        },
        "FS": {
            "tetrahedra": "c++"
        },
        "generate": {
            "hk_from_hr": "python"
        }
    },
    "many_body": {
        "many_body": {
            "triqs": "python",
            "sparse_ir": "julia"
        },
        "self_energy": {
            "sparse_ir": "julia"
        },
        "vertex": {
            "from_susceptibility": "c++"
        },
        "response": {
            "bz_integral": "julia",
            "sparse_ir": "julia"
        }
    },
    "superconductor": {
        "bcs": {
            "lanczos": "c++",  
            "power_iteration": "c++",  
        },
        "eliashberg": {
            "lanczos": "python",  
            "power_iteration": "julia",  
            }
    }
}

if __name__ == "__main__":
    # Run code generators when this file is executed
    print("=" * 60)
    print("FFirefly Category Code Generator")
    print("=" * 60)

    # Import code generation modules
    import sys
    from pathlib import Path

    # Add categories directory to path for imports
    categories_dir = Path(__file__).parent / "categories"
    sys.path.insert(0, str(categories_dir))

    # Import generators
    from write import write
    from write_cmake import update_cmakelists
    # import write_main  # Skip for now - requires CATEGORIES_WITH_TESTS
    # import write_nodes  # Skip for now - not needed for build refactoring

    print("\n1. Generating directory structure and templates...")
    write(CATEGORIES)

    print("\n2. Generating CMakeLists.txt sections...")
    update_cmakelists(CATEGORIES)

    # print("\n3. Generating category node files...")
    # for category in CATEGORIES:
    #     write_nodes.write_node_files(category, CATEGORIES[category])
    # print("   ✓ Generated node.cpp and node.hpp for all categories")

    # print("\n4. Generating main.c sections...")
    # write_main.generate_main_sections(CATEGORIES)
    # print("   ✓ Generated main.c integration code")

    print("\n" + "=" * 60)
    print("Code generation complete!")
    print("=" * 60)
