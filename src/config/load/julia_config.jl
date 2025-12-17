module Config

using PyCall

firefly = pyimport("firefly")
cfg = firefly.config

# Start variable definitions

#[CONTROL]
category::String = cfg.category
calculation::String = cfg.calculation
method::String = cfg.method
outdir::String = cfg.outdir
indir::String = cfg.indir
prefix::String = cfg.prefix
verbosity::String = cfg.verbosity
automatic_file_read::Bool = cfg.automatic_file_read
write_result::Bool = cfg.write_result
filetype::String = cfg.filetype

#[SYSTEM]
interaction::String = cfg.interaction
dimension::Int = cfg.dimension
celltype::String = cfg.celltype
nbnd::Int = cfg.nbnd
nstates::Int = cfg.nstates
fermi_energy::Float64 = cfg.fermi_energy
num_electrons::Float64 = cfg.num_electrons
Temperature::Float64 = cfg.Temperature
onsite_U::Float64 = cfg.onsite_U
cutoff_energy::Float64 = cfg.cutoff_energy
smearing::Float64 = cfg.smearing
mixing::Float64 = cfg.mixing
max_iters::Int = cfg.max_iters
num_solutions::Int = cfg.num_solutions

#[HAMILTONIAN]
hamiltonian::String = cfg.hamiltonian

#[MESH]
k_mesh::Array{Int} = cfg.k_mesh
q_mesh::Array{Int} = cfg.q_mesh
w_pts::Int = cfg.w_pts

#[CELL]
cell::Array{Float64} = cfg.cell

#[BRILLOUIN_ZONE]
brillouin_zone::Array{Float64} = cfg.brillouin_zone

#[BASIS]
states::Vector{String} = cfg.states
positions::Array{Float64} = cfg.positions

#[BANDS]
band::String = cfg.band
eff_mass::Float64 = cfg.eff_mass
t0::Float64 = cfg.t0
t1::Float64 = cfg.t1
t2::Float64 = cfg.t2
t3::Float64 = cfg.t3
t4::Float64 = cfg.t4
t5::Float64 = cfg.t5
t6::Float64 = cfg.t6
t7::Float64 = cfg.t7
t8::Float64 = cfg.t8
t9::Float64 = cfg.t9
t10::Float64 = cfg.t10

#[SUPERCONDUCTOR]
FS_only::Bool = cfg.FS_only
num_eigenvalues_to_save::Int = cfg.num_eigenvalues_to_save
frequency_pts::Int = cfg.frequency_pts
projections::String = cfg.projections

#[RESPONSE]
dynamic::Bool = cfg.dynamic

#[MANY_BODY]
self_consistent::Bool = cfg.self_consistent
# End variable definitions

end
