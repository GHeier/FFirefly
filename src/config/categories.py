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
            "sparse_ir": "julia",
            "triqs": "python",
        },
        "vertex": {
            "from_susceptibility": "c++"
        },
        "response": {
            "tetrahedra": "julia",
            "sparse_ir": "julia"
        },
        "renormalization": {
            "analytic": "c++",
            "from_sigma": "c++"
        }
    },
    "superconductor": {
        "bcs": {
            "convolution": "python",  
            "matrix": "c++",  
        },
        "bcs_w": {
            "hmatrix": "julia",  
        },
        "eliashberg": {
            "convolution": "python",  
            "hmatrix": "julia",  
            "sparse_ir": "julia",  
            }
    }
}

if __name__ == "__main__":
    from write import write_cmake, write_categories
    # Run code generators when this file is executed
    print("=" * 60)
    print("FFirefly Category Code Generator")
    print("=" * 60)

    # Import generators
    # import write_main  # Skip for now - requires CATEGORIES_WITH_TESTS
    # import write_nodes  # Skip for now - not needed for build refactoring

    print("\n1. Generating directory structure and templates...")
    write_categories.write(CATEGORIES)

    print("\n2. Generating CMakeLists.txt sections...")
    write_cmake.update_cmakelists(CATEGORIES)

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
