module Config

using PyCall

firefly = pyimport("firefly")
cfg = firefly.config

# Start variable definitions

#[CONTROL]
category::String = cfg.category
calculation::String = cfg.calculation
outdir::String = cfg.outdir
indir::String = cfg.indir
prefix::String = cfg.prefix
verbosity::String = cfg.verbosity
input_data_file::String = cfg.input_data_file
output_data_file::String = cfg.output_data_file
automatic_file_read::Bool = cfg.automatic_file_read

#[SYSTEM]
interaction::String = cfg.interaction
dimension::Int = cfg.dimension
ibrav::Int = cfg.ibrav
nbnd::Int = cfg.nbnd
natoms::Int = cfg.natoms
fermi_energy::Float64 = cfg.fermi_energy
Temperature::Float64 = cfg.Temperature
onsite_U::Float64 = cfg.onsite_U

#[MESH]
k_mesh::Array{Int} = cfg.k_mesh
q_mesh::Array{Int} = cfg.q_mesh
w_pts::Int = cfg.w_pts

#[CELL]
cell::Array{Float64} = cfg.cell

#[BRILLOUIN_ZONE]
brillouin_zone::Array{Float64} = cfg.brillouin_zone

#[ATOMIC_POSITIONS]
atom::String = cfg.atom
position::Array{Float64} = cfg.position

#[BANDS]
band::Vector{String} = cfg.band

eff_mass::Vector{Float64} = cfg.eff_mass
t0::Vector{Float64} = cfg.t0
t1::Vector{Float64} = cfg.t1
t2::Vector{Float64} = cfg.t2
t3::Vector{Float64} = cfg.t3
t4::Vector{Float64} = cfg.t4
t5::Vector{Float64} = cfg.t5
t6::Vector{Float64} = cfg.t6
t7::Vector{Float64} = cfg.t7
t8::Vector{Float64} = cfg.t8
t9::Vector{Float64} = cfg.t9
t10::Vector{Float64} = cfg.t10

#[SUPERCONDUCTOR]
method::String = cfg.method
FS_only::Bool = cfg.FS_only
bcs_cutoff_frequency::Float64 = cfg.bcs_cutoff_frequency
num_eigenvalues_to_save::Int = cfg.num_eigenvalues_to_save
frequency_pts::Int = cfg.frequency_pts
projections::String = cfg.projections

#[RESPONSE]
dynamic::Bool = cfg.dynamic
# End variable definitions

end
