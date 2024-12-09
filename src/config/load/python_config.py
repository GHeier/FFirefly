import ctypes

lib = ctypes.CDLL('c_config.so')

# Begin Global Variables

#[CONTROL]
category = ''
calculation = ''
outdir = ''
prefix = ''
verbosity = ''
datfile_in = ''
datfile_out = ''

#[SYSTEM]
interaction = ''
dimension = 0
ibrav = 0
nbnd = 0
fermi_energy = 0.0
Temperature = 0.0
onsite_U = 0.0

#[MESH]
k_mesh = [0] * 3
q_mesh = [0] * 3
w_pts = 0

#[CELL]
cell = [[0.0] * 3] * 3
brillouin_zone = [[0.0] * 3] * 3

#[BANDS]
band = [''] * 50

eff_mass = 0.0

#[SUPERCONDUCTOR]
FS_only = False
bcs_cutoff_frequency = 0.0
num_eigenvalues_to_save = 0
frequency_pts = 0

#[RESPONSE]
dynamic = False
# End Global Variables

def load_config():
# Begin Function Declarations

#[CONTROL]
    category_ptr = ctypes.cast(ctypes.addressof(lib.c_category), ctypes.POINTER(ctypes.c_char * 50))
    global category
    category = category_ptr.contents.value.decode('utf-8')
    calculation_ptr = ctypes.cast(ctypes.addressof(lib.c_calculation), ctypes.POINTER(ctypes.c_char * 50))
    global calculation
    calculation = calculation_ptr.contents.value.decode('utf-8')
    outdir_ptr = ctypes.cast(ctypes.addressof(lib.c_outdir), ctypes.POINTER(ctypes.c_char * 50))
    global outdir
    outdir = outdir_ptr.contents.value.decode('utf-8')
    prefix_ptr = ctypes.cast(ctypes.addressof(lib.c_prefix), ctypes.POINTER(ctypes.c_char * 50))
    global prefix
    prefix = prefix_ptr.contents.value.decode('utf-8')
    verbosity_ptr = ctypes.cast(ctypes.addressof(lib.c_verbosity), ctypes.POINTER(ctypes.c_char * 50))
    global verbosity
    verbosity = verbosity_ptr.contents.value.decode('utf-8')
    datfile_in_ptr = ctypes.cast(ctypes.addressof(lib.c_datfile_in), ctypes.POINTER(ctypes.c_char * 50))
    global datfile_in
    datfile_in = datfile_in_ptr.contents.value.decode('utf-8')
    datfile_out_ptr = ctypes.cast(ctypes.addressof(lib.c_datfile_out), ctypes.POINTER(ctypes.c_char * 50))
    global datfile_out
    datfile_out = datfile_out_ptr.contents.value.decode('utf-8')

#[SYSTEM]
    interaction_ptr = ctypes.cast(ctypes.addressof(lib.c_interaction), ctypes.POINTER(ctypes.c_char * 50))
    global interaction
    interaction = interaction_ptr.contents.value.decode('utf-8')
    global dimension
    dimension = ctypes.c_int.in_dll(lib, 'c_dimension').value
    global ibrav
    ibrav = ctypes.c_int.in_dll(lib, 'c_ibrav').value
    global nbnd
    nbnd = ctypes.c_int.in_dll(lib, 'c_nbnd').value
    global fermi_energy
    fermi_energy = ctypes.c_double.in_dll(lib, 'c_fermi_energy').value
    global Temperature
    Temperature = ctypes.c_double.in_dll(lib, 'c_Temperature').value
    global onsite_U
    onsite_U = ctypes.c_double.in_dll(lib, 'c_onsite_U').value

#[MESH]
    global k_mesh
    k_mesh = list((ctypes.c_int * 3).in_dll(lib, 'c_k_mesh'))
    global q_mesh
    q_mesh = list((ctypes.c_int * 3).in_dll(lib, 'c_q_mesh'))
    global w_pts
    w_pts = ctypes.c_int.in_dll(lib, 'c_w_pts').value

#[CELL]
    global cell
    cell = list(((ctypes.c_double * 3) * 3).in_dll(lib, 'c_cell'))
    global brillouin_zone
    brillouin_zone = list(((ctypes.c_double * 3) * 3).in_dll(lib, 'c_brillouin_zone'))

#[BANDS]
    temp = ctypes.POINTER(ctypes.c_char_p)
    global band
    band = temp.in_dll(lib, 'c_band')
    temp = ctypes.c_double * 50
    global eff_mass
    eff_mass = temp.in_dll(lib, 'c_eff_mass').value

#[SUPERCONDUCTOR]
    global FS_only
    FS_only = ctypes.c_bool.in_dll(lib, 'c_FS_only').value
    global bcs_cutoff_frequency
    bcs_cutoff_frequency = ctypes.c_double.in_dll(lib, 'c_bcs_cutoff_frequency').value
    global num_eigenvalues_to_save
    num_eigenvalues_to_save = ctypes.c_int.in_dll(lib, 'c_num_eigenvalues_to_save').value
    global frequency_pts
    frequency_pts = ctypes.c_int.in_dll(lib, 'c_frequency_pts').value

#[RESPONSE]
    global dynamic
    dynamic = ctypes.c_bool.in_dll(lib, 'c_dynamic').value
# End Function Declarations

