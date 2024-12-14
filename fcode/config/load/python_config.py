import ctypes

import os

# Get the absolute path of the current script
current_directory = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(current_directory+'/libc_config.so')

# Begin Function Declarations

#[CONTROL]
global category
category = ctypes.c_char_p.in_dll(lib, 'c_category').value.decode('utf-8')
global calculation
calculation = ctypes.c_char_p.in_dll(lib, 'c_calculation').value.decode('utf-8')
global outdir
outdir = ctypes.c_char_p.in_dll(lib, 'c_outdir').value.decode('utf-8')
global prefix
prefix = ctypes.c_char_p.in_dll(lib, 'c_prefix').value.decode('utf-8')
global verbosity
verbosity = ctypes.c_char_p.in_dll(lib, 'c_verbosity').value.decode('utf-8')
global datfile_in
datfile_in = ctypes.c_char_p.in_dll(lib, 'c_datfile_in').value.decode('utf-8')
global datfile_out
datfile_out = ctypes.c_char_p.in_dll(lib, 'c_datfile_out').value.decode('utf-8')

#[SYSTEM]
global interaction
interaction = ctypes.c_char_p.in_dll(lib, 'c_interaction').value.decode('utf-8')
global dimension
dimension = ctypes.c_int.in_dll(lib, 'c_dimension').value
global ibrav
ibrav = ctypes.c_int.in_dll(lib, 'c_ibrav').value
global nbnd
nbnd = ctypes.c_int.in_dll(lib, 'c_nbnd').value
global fermi_energy
fermi_energy = ctypes.c_float.in_dll(lib, 'c_fermi_energy').value
global Temperature
Temperature = ctypes.c_float.in_dll(lib, 'c_Temperature').value
global onsite_U
onsite_U = ctypes.c_float.in_dll(lib, 'c_onsite_U').value

#[MESH]
global k_mesh
k_mesh = list((ctypes.c_int * 3).in_dll(lib, 'c_k_mesh'))
global q_mesh
q_mesh = list((ctypes.c_int * 3).in_dll(lib, 'c_q_mesh'))
global w_pts
w_pts = ctypes.c_int.in_dll(lib, 'c_w_pts').value

#[CELL]
global cell
cell = [[(((ctypes.c_float * 3) * 3).in_dll(lib, 'c_cell'))[i][j] for j in range(3)] for i in range(3)]
global brillouin_zone
brillouin_zone = [[(((ctypes.c_float * 3) * 3).in_dll(lib, 'c_brillouin_zone'))[i][j] for j in range(3)] for i in range(3)]

#[BANDS]
global band
band = ctypes.c_char_p.in_dll(lib, 'c_band').value.decode('utf-8')
global eff_mass
eff_mass = ctypes.c_float.in_dll(lib, 'c_eff_mass').value
global t0
t0 = ctypes.c_float.in_dll(lib, 'c_t0').value
global t1
t1 = ctypes.c_float.in_dll(lib, 'c_t1').value
global t2
t2 = ctypes.c_float.in_dll(lib, 'c_t2').value
global t3
t3 = ctypes.c_float.in_dll(lib, 'c_t3').value
global t4
t4 = ctypes.c_float.in_dll(lib, 'c_t4').value
global t5
t5 = ctypes.c_float.in_dll(lib, 'c_t5').value
global t6
t6 = ctypes.c_float.in_dll(lib, 'c_t6').value
global t7
t7 = ctypes.c_float.in_dll(lib, 'c_t7').value
global t8
t8 = ctypes.c_float.in_dll(lib, 'c_t8').value
global t9
t9 = ctypes.c_float.in_dll(lib, 'c_t9').value
global t10
t10 = ctypes.c_float.in_dll(lib, 'c_t10').value

#[SUPERCONDUCTOR]
global FS_only
FS_only = ctypes.c_bool.in_dll(lib, 'c_FS_only').value
global bcs_cutoff_frequency
bcs_cutoff_frequency = ctypes.c_float.in_dll(lib, 'c_bcs_cutoff_frequency').value
global num_eigenvalues_to_save
num_eigenvalues_to_save = ctypes.c_int.in_dll(lib, 'c_num_eigenvalues_to_save').value
global frequency_pts
frequency_pts = ctypes.c_int.in_dll(lib, 'c_frequency_pts').value

#[RESPONSE]
global dynamic
dynamic = ctypes.c_bool.in_dll(lib, 'c_dynamic').value
# End Function Declarations

