# The format is as follows:
# SECTION_NAME: {
#     VARIABLE_NAME: DEFAULT_VALUE,
#     ...
# }
# Follow this format when adding new sections and variables.
# Default values are the values used when the variable is not defined in the input file
# Cannot pass strings and bools through arrays

import numpy as np

ALL = {
    "CONTROL": {
        "category": "test",
        "calculation": "test",
        "method": "none",
        "outdir": "./",
        "debug": False,
        "prefix": "sample",
        "verbosity": "low",
        "automatic_file_read": True,
        "write_result": True,
        "filetype": 'h5',
    },
    "SYSTEM": {
        "interaction": "none",
        "dimension": 3,
        "celltype": "",
        "nbnd": 0,
        "fermi_energy": 0.0,
        "num_electrons": 0.0,
        "mu_from_n": False,
        "Temperature": 0.0,
        "cutoff_energy": 0.05,
        "smearing": 0.02,
        "mixing": 0.02,
        "max_iters": 100,
        "qp_weight": 1.0,
    },
    "HAMILTONIAN": {"hamiltonian": "tight_binding"},
    "HUBBARD": {"U0": 0.0, "U1": 0.0, "J0": 0.0, "J1": 0.0},
    "MESH": {"k_mesh": [10, 10, 10], "q_mesh": [10, 10, 10], "w_pts": 100},
    "CELL": {"cell": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]},
    "BRILLOUIN_ZONE": {
        "brillouin_zone": [
            [2 * np.pi, 0.0, 0.0],
            [0.0, 2 * np.pi, 0.0],
            [0.0, 0.0, 2 * np.pi],
        ]
    },
    "BASIS": {"states": ["H"], "positions": [[0.0, 0.0, 0.0]]},  # Max 50 states
    "BANDS": {
        "band": "fermi_gas",
        "eff_mass": 1.0,
        "t0": 1.0,
        "t1": 0.0,
        "t2": 0.0,
        "t3": 0.0,
        "t4": 0.0,
        "t5": 0.0,
        "t6": 0.0,
        "t7": 0.0,
        "t8": 0.0,
        "t9": 0.0,
        "t10": 0.0,
    },
    "SUPERCONDUCTOR": {
        "FS_only": True,
        "num_eigenvalues_to_save": 5,
        "frequency_pts": 0,
        "projections": "",
    },
    "RESPONSE": {"dynamic": False},
    "MANY_BODY": {
        "self_consistent": False,
        "impurity_solver": "IPT"
    },
}


from write import write_c, write_cpp, write_f90, write_py, write_julia

write_c.write_c(ALL)
write_c.write_c_header(ALL)
write_cpp.write_cpp(ALL)
write_cpp.write_cpp_header(ALL)
write_f90.write_f90(ALL)
write_py.write_py(ALL)
write_julia.write_julia(ALL)
