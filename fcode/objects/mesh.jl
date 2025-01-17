module IRMesh

using FFTW
using LinearAlgebra
using Roots
using SparseIR
import SparseIR: Statistics, value, valueim

export Mesh
export tau_to_wn, wn_to_tau, k_to_r, r_to_k

"""
Holding struct for k-mesh and sparsely sampled imaginary time 'tau' / Matsubara frequency 'iw_n' grids.
Additionally we defines the Fourier transform routines 'r <-> k'  and 'tau <-> l <-> wn'.
 """
struct Mesh
    nk1         ::Int64
    nk2         ::Int64
    nk          ::Int64
    ek          ::Array{Float64,2}
    iw0_f       ::Int64
    iw0_b       ::Int64
    fnw         ::Int64
    fntau       ::Int64
    bnw         ::Int64
    bntau       ::Int64
    IR_basis_set::FiniteTempBasisSet
end

"""Initiarize function"""
function Mesh(
        nk1         ::Int64,
        nk2         ::Int64,
        IR_basis_set::FiniteTempBasisSet,
        )::Mesh

    nk::Int64 = nk1*nk2

    # Compute Hamiltonian
    ek = Array{ComplexF64,2}(undef, nk1, nk2)
    for iy in 1:nk2, ix in 1:nk1
        kx::Float64 = (2*π*(ix-1))/nk1
        ky::Float64 = (2*π*(iy-1))/nk2
        ek[ix, iy] = -2.0*(cos(kx)+cos(ky))
    end

    # lowest Matsubara frequency index
    iw0_f = findall(x->x==FermionicFreq(1), IR_basis_set.smpl_wn_f.sampling_points)[1]
    iw0_b = findall(x->x==BosonicFreq(0), IR_basis_set.smpl_wn_b.sampling_points)[1]

    # the number of sampling point for fermion and boson
    fnw   = length(IR_basis_set.smpl_wn_f.sampling_points)
    fntau = length(IR_basis_set.smpl_tau_f.sampling_points)
    bnw   = length(IR_basis_set.smpl_wn_b.sampling_points)
    bntau = length(IR_basis_set.smpl_tau_b.sampling_points)

    # Return
    Mesh(nk1, nk2, nk, ek, iw0_f, iw0_b, fnw, fntau, bnw, bntau, IR_basis_set)
end

function smpl_obj(mesh::Mesh, statistics::SparseIR.Statistics)
    """ Return sampling object for given statistic """
    if statistics == Fermionic()
        smpl_tau = mesh.IR_basis_set.smpl_tau_f
        smpl_wn  = mesh.IR_basis_set.smpl_wn_f
    elseif statistics == Bosonic()
        smpl_tau = mesh.IR_basis_set.smpl_tau_b
        smpl_wn  = mesh.IR_basis_set.smpl_wn_b
    end
    return smpl_tau, smpl_wn
end

"""Fourier transformation"""
function tau_to_wn(mesh::Mesh, statistics::T, obj_tau) where {T <: SparseIR.Statistics}
    """ Fourier transform from tau to iw_n via IR basis """
    smpl_tau, smpl_wn = smpl_obj(mesh, statistics)

    obj_l = fit(smpl_tau, obj_tau, dim=1)
    obj_wn = evaluate(smpl_wn, obj_l, dim=1)
    return obj_wn
end

function wn_to_tau(mesh::Mesh, statistics::Statistics, obj_wn)
    """ Fourier transform from iw_n to tau via IR basis """
    smpl_tau, smpl_wn = smpl_obj(mesh, statistics)

    obj_l   = fit(smpl_wn, obj_wn, dim=1)
    obj_tau = evaluate(smpl_tau, obj_l, dim=1)
    return obj_tau
end

function k_to_r(mesh::Mesh, obj_k)
    """ Fourier transform from k-space to real space """
    obj_r = fft(obj_k,[2,3])
    return obj_r
end

function r_to_k(mesh::Mesh, obj_r)
    """ Fourier transform from real space to k-space """
    obj_k = ifft(obj_r,[2,3])/mesh.nk
    return obj_k
end

end
