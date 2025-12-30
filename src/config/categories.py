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
    import os
    import sys
    import importlib.util
    from categories.write import write
    write(CATEGORIES)
