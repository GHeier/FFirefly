module Linearized_Eliashberg2
using LaTeXStrings

using FFTW
using LinearAlgebra
using Roots
using SparseIR
import SparseIR: Statistics, value, valueim
include("../objects/mesh.jl")
using .IRMesh

using Firefly
cfg = Firefly.Config

### System parameters
T    = cfg.Temperature
beta = 1/T    # inverse temperature
U    = cfg.onsite_U
mu   = cfg.fermi_energy
BZ   = cfg.brillouin_zone
wc   = cfg.cutoff_energy
outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype

### Numerical parameters
nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim == 2
    nz = 1
end
nw = cfg.w_pts
nk        = nx*ny*nz
sfc_tol   = 1e-4      # accuracy for self-consistent iteration
maxiter   = 30        # maximal number of iterations in self-consistent cycle

bcs_debug = false
;

"""
Solver struct for solving the linearized gap equation using the power method.
It takes FLEX results as an input.
"""
mutable struct LinearizedGapSolver
    mesh      ::Mesh
    gkio      ::Array{ComplexF64,4}
    V_singlet ::Array{ComplexF64,4}
    delta     ::Array{ComplexF64,4}
    frit      ::Array{ComplexF64,4}
    iw        ::Vector{ComplexF64}
    maxiter   ::Int64
    sfc_tol   ::Float64
    verbose   ::Bool
    lam       ::Float64
end

function LinearizedGapSolver(
        mesh       ::Mesh,
        gkio       ::Array{ComplexF64,4};
        iw         ::Vector{ComplexF64},
        V          ::Array{ComplexF64,4},
        maxiter    ::Int64  =50,
        sfc_tol    ::Float64=1e-4,
        verbose    ::Bool   =true
        )::LinearizedGapSolver

    ## Initialize necessary quantities from FLEX loop

    maxiter = maxiter
    sfc_tol = sfc_tol
    verbose = verbose

    ## Initialize trial gap function
    # Here we focus on a d-wave symmetric solution
    delta = Array{ComplexF64}(undef, mesh.fnw, nx, ny, nz)
    frit  = Array{ComplexF64}(undef, mesh.fntau, nx, ny, nz)
    for iy in 1:ny, ix in 1:nx, iw in 1:mesh.fnw
        kx::Float64 = 2*π* ((ix-1)/nx - 0.5)
        ky::Float64 = 2*π* ((iy-1)/ny - 0.5)
        delta[iw,ix,iy] = cos(kx) - cos(ky)
    end
    if bcs_debug
        delta = 1.0 * ones(ComplexF64, nw, nx, ny, nz)
        frit = Array{ComplexF64}(undef, nw, nx, ny, nz)
        zero_out_beyond_wc_iw!(delta, iw)
    end

    #normalize initial guess
    normalize!(delta)

    # Initialize interaction
    V_singlet = V

    ## Initialize eigenvalue
    lam::Float64 = 0.0
    gap_solver = LinearizedGapSolver(mesh, gkio, V_singlet, delta, frit, iw, maxiter, sfc_tol, verbose, lam)
end

function solve(gap_solver::LinearizedGapSolver)
    """ Solving instance to find eigenvalue from power method """
    for it in 1:gap_solver.maxiter
        lam_old = gap_solver.lam
        delta_old = copy(gap_solver.delta)

        frit_calc(gap_solver)
        deltarit = gap_solver.V_singlet .* gap_solver.frit

        if bcs_debug
            gap_solver.delta .= fftshift(ifft(deltarit)) / (nk * beta)
            zero_out_beyond_wc_iw!(gap_solver.delta, gap_solver.iw)
        else
            gap_solver.delta .= rtau_to_kw(deltarit, 'F', gap_solver.mesh)
        end

        # calculate eigenvalue
        gap_solver.lam = sum(real.(conj.(gap_solver.delta).* delta_old))
        normalize!(gap_solver.delta)
        #deflate_v!(gap_solver.delta, oldvecs)

        if gap_solver.verbose
            println(it, '\t', gap_solver.lam, '\t', abs(gap_solver.lam-lam_old))
        end
        if abs(gap_solver.lam-lam_old) < gap_solver.sfc_tol
            break
        end
    end
end


function frit_calc(gap_solver::LinearizedGapSolver)
    """ Calculate (linearized) anomalous Green function F = |G|^2 * delta for evaluating the gap equation """
    fkio = - gap_solver.gkio.*conj(gap_solver.gkio).*gap_solver.delta
    if bcs_debug
        gap_solver.frit = fft(fkio)
    else
        gap_solver.frit = kw_to_rtau(fkio, 'F', gap_solver.mesh)
    end
end

function deflate_v!(v, eigvecs)
    for phi_i in eigvecs
        proj = (dot(vec(phi_i), vec(v)) / dot(vec(phi_i), vec(phi_i))) * phi_i
        v .-= proj
    end
    nv = norm(vec(v))
    if nv > 0
        v ./= nv
    else
        error("Deflation resulted in zero vector")
    end
end

function create_energy_mesh(band, with_sigma=true, iw = 0, Sigma = 0)
    nw = length(iw)
    e_arr = zeros(ComplexF64, nw, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nw
        kvec = get_kvec(i, j, k)
        if with_sigma
            e_arr[l, i, j, k] = band(1, kvec) - mu + Sigma(kvec, imag(iw[l]))
        else
            e_arr[l, i, j, k] = band(1, kvec) - mu
        end
    end
    return e_arr
end

function create_G(band, iw, Sigma)
    nw = length(iw)
    e_arr = zeros(ComplexF64, nw, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nw
        kvec = get_kvec(i, j, k)
        if bcs_debug
            e_arr[l, i, j, k] = band(1, kvec) - mu 
        else
        e_arr[l, i, j, k] = band(1, kvec) - mu + Sigma(kvec, imag(iw[l]))
        end
    end
    return (1 ./ (iw .- e_arr))
end

function fill_V_rt(vertex, iv, mesh)
    bnw = length(iv)
    V = Array{ComplexF64}(undef, bnw, nx, ny, nz)
    if bcs_debug
        V = -1.0 * ones(ComplexF64, nw, nx, ny, nz)
        return fft(V)
    end
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:bnw
        q = get_kvec(i, j, k)
        w = imag(iv[l])
        #V[l, i, j, k] = 0.5 * (vertex(q, w) + vertex(-q, w))
        V[l, i, j, k] = vertex(q, w)
    end
    return kw_to_rtau(V, 'B', mesh)
end

function eigenvalue_computation()
    # initialize calculation
    band = Bands()
    e_arr = create_energy_mesh(band, false)

    mesh = IR_Mesh(real.(e_arr))
    iw, iv = get_iw_iv(mesh)
    if bcs_debug
        iw, iv = get_iw_iv_bcs()
        N = Firefly.Field_R(outdir * prefix * "_DOS." * filetype)
        initial_expected = log(1.134 * wc * beta) * N(mu)
        println("Initial Expected eig: $initial_expected")
    end
    iw = ComplexF64.(iw)
    iv = ComplexF64.(iv)

    self_energy = Self_Energy()
    vertex = Vertex()
    gkio = create_G(band, iw, self_energy)
    V_rt = fill_V_rt(vertex, iv, mesh)

    gap_solver = LinearizedGapSolver(mesh, gkio; iw=iw, V=V_rt, maxiter=maxiter, sfc_tol=sfc_tol)
    solve(gap_solver)
    println("The superconducting eigenvalue is lambda_d=",gap_solver.lam)
end

function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
end

function zero_out_beyond_wc_iw!(A, iw)
    if size(A, 1) != length(iw)
        println("Zero_Out: Incorrect array mapping")
    end
    for i in 1:length(iw)
        if abs(imag(iw[i])) > wc 
            A[i, :, :, :] .= 0
        end
    end
end

function get_iw_iv_bcs()
    iw = Array{ComplexF32}(undef, nw)
    iv = Array{ComplexF32}(undef, nw)
    for i in 1:nw
        n = -nw / 2 + i
        iw[i] = pi / beta * (2 * n + 1) * im
    end
    for i in 1:nw
        n = -nw / 2 + i
        iv[i] = pi / beta * (2 * n) * im
    end
    println("Fermionic frequency from ", minimum(imag.(iw)), " to ", maximum(imag.(iw)))
    println("Bosonic frequency from ", minimum(imag.(iv)), " to ", maximum(imag.(iv)))
    return iw, iv 
end
end # module
