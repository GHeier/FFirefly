module Config

# NOTE: This module provides access to config variables.
# PyCall is NOT used to avoid conflicts with the Python interpreter
# that is already initialized by the C code in start_python().
# Instead, values are loaded from default or through C bindings.

# Start variable definitions

#[CONTROL]
category::String = "test"
calculation::String = "test"
method::String = "none"
outdir::String = "./"
indir::String = "./"
prefix::String = "sample"
verbosity::String = "low"
automatic_file_read::Bool = true
write_result::Bool = true
filetype::String = "h5"

#[SYSTEM]
interaction::String = "none"
dimension::Int = 3
celltype::String = ""
nbnd::Int = 0
nstates::Int = 0
fermi_energy::Float64 = 0.0
num_electrons::Float64 = 0.0
Temperature::Float64 = 0.0
onsite_U::Float64 = 0.0
cutoff_energy::Float64 = 0.05
smearing::Float64 = 0.02
mixing::Float64 = 0.02
max_iters::Int = 100
num_solutions::Int = 5

#[HAMILTONIAN]
hamiltonian::String = "tight_binding"

#[MESH]
k_mesh::Array{Int} = [10, 10, 10]
q_mesh::Array{Int} = [10, 10, 10]
w_pts::Int = 100

#[CELL]
cell::Vector{Vector{Float64}} = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

#[BRILLOUIN_ZONE]
brillouin_zone::Vector{Vector{Float64}} = [[6.283185307179586, 0.0, 0.0], [0.0, 6.283185307179586, 0.0], [0.0, 0.0, 6.283185307179586]]

#[BASIS]
states::Vector{String} = ["H"]
positions::Vector{Vector{Float64}} = [[0.0, 0.0, 0.0]]

#[BANDS]
band::String = "fermi_gas"
eff_mass::Float64 = 1.0
t0::Float64 = 1.0
t1::Float64 = 0.0
t2::Float64 = 0.0
t3::Float64 = 0.0
t4::Float64 = 0.0
t5::Float64 = 0.0
t6::Float64 = 0.0
t7::Float64 = 0.0
t8::Float64 = 0.0
t9::Float64 = 0.0
t10::Float64 = 0.0

#[SUPERCONDUCTOR]
FS_only::Bool = true
num_eigenvalues_to_save::Int = 0
frequency_pts::Int = 0
projections::String = ""

#[RESPONSE]
dynamic::Bool = false

#[MANY_BODY]
self_consistent::Bool = false
# End variable definitions

end
