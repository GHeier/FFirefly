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
        "outdir": "./",
        "indir": "./",
        "prefix": "sample",
        "verbosity": "low",
        "input_data_file": "input.dat",
        "output_data_file": "output.dat",
        "automatic_file_read": True,
    },
    "SYSTEM": {
        "interaction": "none",
        "dimension": 3,
        "celltype": "",
        "nbnd": 0,
        "natoms": 0,
        "fermi_energy": 0.0,
        "Temperature": 0.0,
        "onsite_U": 0.0,
    },
    "MESH": {"k_mesh": [10, 10, 10], "q_mesh": [10, 10, 10], "w_pts": 100},
    "CELL": {"cell": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]},
    "BRILLOUIN_ZONE": {
        "brillouin_zone": [
            [2 * np.pi, 0.0, 0.0],
            [0.0, 2 * np.pi, 0.0],
            [0.0, 0.0, 2 * np.pi],
        ]
    },
    "ATOMIC_POSITIONS": {"atom": "X", "position": [0.0, 0.0, 0.0]},
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
        "method": "none",
        "FS_only": True,
        "bcs_cutoff_frequency": 0.05,
        "num_eigenvalues_to_save": 1,
        "frequency_pts": 5,
        "projections": "",
    },
    "RESPONSE": {"dynamic": False},
}


from write import write_c, write_cpp, write_f90, write_py, write_julia

write_c.write_c(ALL)
write_c.write_c_header(ALL)
write_cpp.write_cpp(ALL)
write_cpp.write_cpp_header(ALL)
write_f90.write_f90(ALL)
write_py.write_py(ALL)
write_julia.write_julia(ALL)
