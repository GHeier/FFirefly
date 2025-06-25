module IRMesh

using FFTW
using LinearAlgebra
using Roots
using SparseIR
import SparseIR: Statistics, value, valueim
using PyCall

using Firefly

cfg = Firefly.Config

BZ = cfg.brillouin_zone
dim = cfg.dimension
beta = 1 / cfg.Temperature

export Mesh, get_iw_iv
export tau_to_wn, wn_to_tau, k_to_r, r_to_k, kw_to_rtau, rtau_to_kw, IR_Mesh
export k_to_r!, r_to_k!, kw_to_rtau!, rtau_to_kw!, IR_Mesh

"""
Holding struct for k-mesh and sparsely sampled imaginary time 'tau' / Matsubara frequency 'iw_n' grids.
Additionally we defines the Fourier transform routines 'r <-> k'  and 'tau <-> l <-> wn'.
 """
struct Mesh
    nk1         ::Int64
    nk2         ::Int64
    nk          ::Int64
    iw0_f       ::Int64
    iw0_b       ::Int64
    fnw         ::Int64
    fntau       ::Int64
    bnw         ::Int64
    bntau       ::Int64
    IR_basis_set::FiniteTempBasisSet
    dimension   ::Int64
end

"""Initialize function"""
function Mesh(
        IR_basis_set::FiniteTempBasisSet,
        nk1         ::Int64 = 0,
        nk2         ::Int64 = 0,
        nk3         ::Int64 = 0,
        )::Mesh

    nk::Int64 = nk1 * nk2 * nk3
    dimension = 3
    if nk3 == 1
        dimension = 2
    end
    if nk2 == 1
        dimension = 1
    end
    if nk1 == 1
        dimension = 0
    end
    println("IRMesh: nk = ", nk, " dimension = ", dimension)

    # lowest Matsubara frequency index
    iw0_f = findall(x->x==FermionicFreq(1), IR_basis_set.smpl_wn_f.sampling_points)[1]
    iw0_b = findall(x->x==BosonicFreq(0), IR_basis_set.smpl_wn_b.sampling_points)[1]

    # the number of sampling point for fermion and boson
    fnw   = length(IR_basis_set.smpl_wn_f.sampling_points)
    fntau = length(IR_basis_set.smpl_tau_f.sampling_points)
    bnw   = length(IR_basis_set.smpl_wn_b.sampling_points)
    bntau = length(IR_basis_set.smpl_tau_b.sampling_points)

    # Return
    Mesh(nk1, nk2, nk, iw0_f, iw0_b, fnw, fntau, bnw, bntau, IR_basis_set, dimension)
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
    #max_l = maximum(abs.(obj_l))
    #for i in 1:length(obj_l)
    #    if abs(obj_l[i]) < 1e-5 * max_l
    #        obj_l[i] = 0
    #    end
    #end
    obj_tau = evaluate(smpl_tau, obj_l, dim=1)
    return obj_tau
end

function tau_to_wn!(obj_wn, statistics::T, mesh::Mesh, obj_tau) where {T <: SparseIR.Statistics}
    """ Fourier transform from tau to iw_n via IR basis """
    smpl_tau, smpl_wn = smpl_obj(mesh, statistics)

    obj_l = fit(smpl_tau, obj_tau, dim=1)
    obj_wn .= evaluate(smpl_wn, obj_l, dim=1)
end

function wn_to_tau!(obj_tau, statistics::Statistics, mesh::Mesh, obj_wn)
    smpl_tau, smpl_wn = smpl_obj(mesh, statistics)

    obj_l = fit(smpl_wn, obj_wn, dim=1)  # This allocates and must be temporary

    # Truncation (in-place) on obj_l
    max_l = maximum(abs.(obj_l))
    for i in eachindex(obj_l)
        if abs(obj_l[i]) < 1e-5 * max_l
            obj_l[i] = 0
        end
    end

    # Evaluate directly into preallocated output array
    obj_tau .= evaluate(smpl_tau, obj_l, dim=1)
end

function tau_to_wn_c(mesh::Mesh, statistics::Statistics, obj_tau)
    """ Custom Fourier transform from tau to iw_n """
    obj_wn = ifft(obj_tau, [1])
    return obj_wn 
end

function wn_to_tau_c(mesh::Mesh, statistics::Statistics, obj_wn)
    """ Custom Fourier transform from iw_n to tau """
    obj_tau = fft(obj_wn, [1])
    return obj_tau
end

function k_to_r(obj_k, mesh)
    """ Fourier transform from k-space to real space """
    if mesh.dimension == 1
        obj_r = fft(obj_k,[2])
    elseif mesh.dimension == 2
        obj_r = fft(obj_k,[2,3])
    elseif mesh.dimension == 3
        obj_r = fft(obj_k,[2,3,4])
    end
    return obj_r 
end

function r_to_k(obj_r, mesh)
    """ Fourier transform from real space to k-space """
    if mesh.dimension == 1
        obj_k = ifft(obj_r,[2])
    elseif mesh.dimension == 2
        obj_k = ifft(obj_r,[2,3])
    elseif mesh.dimension == 3
        obj_k = ifft(obj_r,[2,3,4])
    end
    return obj_k / mesh.nk
end

function k_to_r!(obj_r::AbstractArray, obj_k::AbstractArray, mesh::Mesh)
    obj_r .= fft(obj_k, ntuple(i -> i+1, mesh.dimension))
end

function r_to_k!(obj_k, obj_r, mesh)
    obj_k .= ifft(obj_r, ntuple(i -> i+1, mesh.dimension)) ./ mesh.nk
end

function kw_to_rtau!(out, obj_k, particle_type::Char, mesh)
    if particle_type == 'F'
        wn_to_tau!(out, Fermionic(), mesh, obj_k)  # Also needs to be in-place
    elseif particle_type == 'B'
        wn_to_tau!(out, Bosonic(), mesh, obj_k)
    end
    k_to_r!(out, out, mesh)  # This must be an in-place version
end

function rtau_to_kw!(out, obj_r, particle_type, mesh)
    if particle_type == 'F'
        tau_to_wn!(out, Fermionic(), mesh, obj_r)  # Also needs to be in-place
    elseif particle_type == 'B'
        tau_to_wn!(out, Bosonic(), mesh, obj_r)  # Also needs to be in-place
    end
    r_to_k!(out, out, mesh)
end

function kw_to_rtau(obj_k, particle_type, mesh)
    temp = k_to_r(obj_k, mesh)
    if particle_type == 'F'
        return wn_to_tau(mesh, Fermionic(), temp)
    elseif particle_type == 'B'
        return wn_to_tau(mesh, Bosonic(), temp)
    end
end

function rtau_to_kw(obj_r, particle_type, mesh)
    temp = r_to_k(obj_r, mesh)
    if particle_type == 'F'
        return tau_to_wn(mesh, Fermionic(), temp)
    elseif particle_type == 'B'
        return tau_to_wn(mesh, Bosonic(), temp)
    end
end

function get_bandwidth()
    maxval = -1000
    minval = 1000
    nz = 200
    if dim == 2
        nz = 1
    end
    for i in 1:200, j in 1:200, k in 1:nz
        kvec = BZ * [i / 200, j / 200, k / 200]
        eps = epsilon(1, kvec)
        maxval = max(maxval, eps)
        minval = min(minval, eps)
    end
    return maxval - minval
end

function IR_Mesh(energy_mesh::Array{Float32, 4}, IR_tol = 1e-10)
    beta = 1 / cfg.Temperature
    _, nx, ny, nz = size(energy_mesh)
    wmax = maximum(energy_mesh) - minimum(energy_mesh)
    IR_basis_set = FiniteTempBasisSet(beta, Float64(1.2*wmax), IR_tol)
    return Mesh(IR_basis_set, nx, ny, nz)
    if dim == 1
        mesh = Mesh(IR_basis_set, nx, 0, 0)
    elseif dim == 2
        mesh = Mesh(IR_basis_set, nx, ny, 0)
    elseif dim == 3
        mesh = Mesh(IR_basis_set, nx, ny, nz)
    end
    println("mesh dimension is ", mesh.dimension)
    return mesh
end


function get_iw_iv(mesh)
    fnw, bnw, fntau, bntau = mesh.fnw, mesh.bnw, mesh.fntau, mesh.bntau
    println("fnw: ", fnw, " bnw: ", bnw, " fntau: ", fntau, " bntau: ", bntau)
    println("Created IRMesh")
    iw = Array{ComplexF32}(undef, fnw)
    iv = Array{ComplexF32}(undef, bnw)
    for i in 1:fnw
        iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
    end
    for i in 1:bnw
        iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
    end
    println("Fermionic frequency from ", minimum(imag.(iw)), " to ", maximum(imag.(iw)))
    println("Bosonic frequency from ", minimum(imag.(iv)), " to ", maximum(imag.(iv)))
    return iw, iv 
end


end
