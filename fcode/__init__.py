from .fmodule import *

#from .algorithms import sparse_ir_mesh
from .config.load import config
from .sample import something

config.load_config()




placeholder = Config()
# Begin python->c++ interface

#[CONTROL]
placeholder.category = config.category
placeholder.calculation = config.calculation
placeholder.outdir = config.outdir
placeholder.prefix = config.prefix
placeholder.verbosity = config.verbosity
placeholder.datfile_in = config.datfile_in
placeholder.datfile_out = config.datfile_out

#[SYSTEM]
placeholder.interaction = config.interaction
placeholder.dimension = config.dimension
placeholder.ibrav = config.ibrav
placeholder.nbnd = config.nbnd
placeholder.fermi_energy = config.fermi_energy
placeholder.Temperature = config.Temperature
placeholder.onsite_U = config.onsite_U

#[MESH]
placeholder.k_mesh = config.k_mesh
placeholder.q_mesh = config.q_mesh
placeholder.w_pts = config.w_pts

#[CELL]
placeholder.cell = config.cell

#[BRILLOUIN_ZONE]
placeholder.brillouin_zone = config.brillouin_zone

#[BANDS]
placeholder.band = config.band
placeholder.eff_mass = config.eff_mass
placeholder.t0 = config.t0
placeholder.t1 = config.t1
placeholder.t2 = config.t2
placeholder.t3 = config.t3
placeholder.t4 = config.t4
placeholder.t5 = config.t5
placeholder.t6 = config.t6
placeholder.t7 = config.t7
placeholder.t8 = config.t8
placeholder.t9 = config.t9
placeholder.t10 = config.t10

#[SUPERCONDUCTOR]
placeholder.method = config.method
placeholder.FS_only = config.FS_only
placeholder.bcs_cutoff_frequency = config.bcs_cutoff_frequency
placeholder.num_eigenvalues_to_save = config.num_eigenvalues_to_save
placeholder.frequency_pts = config.frequency_pts

#[RESPONSE]
placeholder.dynamic = config.dynamic
# End python->c++ interface

py_to_cpp_config(placeholder)
