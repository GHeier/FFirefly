# The format is as follows:
# SECTION_NAME: {
#     VARIABLE_NAME: DEFAULT_VALUE,
#     ...
# }
# Follow this format when adding new sections and variables.
# Default values are the values used when the variable is not defined in the input file
# Cannot pass strings and bools through arrays

import os
import sys
import numpy as np

ALL = {
    'CONTROL': {
        'category': 'test',
        'calculation': 'test',
        'outdir': './output',
        'prefix': 'sample',
        'verbosity': 'high',
        'datfile_in': 'input.dat',
        'datfile_out': 'output.dat'
    },
    'SYSTEM': {
        'interaction': 'none',
        'dimension': 3,
        'ibrav': 0,
        'fermi_energy': 0.0,
        'onsite_U': 0.0
    },
    'MESH': {
        'k_mesh': [10, 10, 10],
        'q_mesh': [10, 10, 10],
        'w_pts': 100
    },
    'CELL': {
        'cell': [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        'brillouin_zone': [[2 * np.pi, 0.0, 0.0], [0.0, 2 * np.pi, 0.0], [0.0, 0.0, 2 * np.pi]]
    },
    'BANDS': {
        'bands': 'fermi_gas',
        'eff_mass': 1.0
    },
    'SUPERCONDUCTOR': {
        'FS_only': True,
        'bcs_cutoff_frequency': 0.05,
        'num_eigenvalues_to_save': 1,
        'frequency_pts': 5
    },
    'RESPONSE': {
        'dynamic': False
    }
}


sys.path.append(os.path.join(os.path.dirname(__file__), 'write'))
from write_c import *
from write_cpp import *
from write_f90 import *
write_c(ALL)
write_c_header(ALL)
write_cpp(ALL)
write_cpp_header(ALL)
write_f90(ALL)
