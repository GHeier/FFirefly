module ManyBodyLoop

using Random
using PencilFFTs
using MPI
using LoopVectorization
using FFTW
using LinearAlgebra
using Roots
using SparseIR
import SparseIR: Statistics, value, valueim
using Printf
using Interpolations
include("../objects/mesh.jl")
using .IRMesh
using Firefly
cfg = Firefly.Config

outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype

scf = cfg.self_consistent 
scf_tol = 1e-4
mix = 0.2

nx, ny, nz = cfg.k_mesh
nqx, nqy, nqz = cfg.q_mesh
dim = cfg.dimension
if dim == 2
    nz = 1
    nqz = 1
end
nk = nx * ny * nz
nbnd = cfg.nbnd
mu = cfg.fermi_energy
U = cfg.onsite_U
BZ = cfg.brillouin_zone
beta = 1 / cfg.Temperature


### System parameters
T    = cfg.Temperature
beta = 1/T    # inverse temperature
#n    = 0.85   # electron filling, here per spin per lattice site (n=1: half filling)
nbnd = cfg.nbnd
mu = cfg.fermi_energy
U = cfg.onsite_U
BZ = cfg.brillouin_zone

### Numerical parameters
nk1, nk2, nk3  = cfg.k_mesh
dim = cfg.dimension
if dim == 2 nk3 = 1 end
nk        = nk1*nk2*nk3
sfc_tol   = 1e-4      # accuracy for self-consistent iteration
maxiter   = 30        # maximal number of iterations in self-consistent cycle
if scf == false maxiter = 1 end
mix       = 0.2       # mixing parameter for new 
U_maxiter = 50       # maximal number of iteration steps in U renormalization loop

interaction = cfg.interaction

function fill_energy_mesh_mpi!(band, ek)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i - 1, j - 1, k - 1, nx, ny, nz)
        ek[1, i, j, k] = band(1, kvec)
    end
    return ek
end

function fill_energy_mesh(band)
    ek = Array{Float32}(undef, 1, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i - 1, j - 1, k - 1, nx, ny, nz)
        ek[1, i, j, k] = band(1, kvec)
    end
    return ek
end

function get_kvec(ix, iy, iz, nx, ny, nz)
    kvec = [ix / (nx) - 0.0, iy / (ny) - 0.0, iz / (nz) - 0.0] 
    if dim < 3
        kvec[3] = 0.0
    elseif dim < 2
        kvec[2] = 0.0
    end
    kvec = BZ * kvec
    return kvec
end

"""
Solver struct to calculate the FLEX loop self-consistently.
After initializing the Solver by `solver = ManyBodySolver(mesh, beta, U, n, sigma_init, sfc_tol, maxiter, U_maxiter, mix)'
it can be run by `solve(solver)`.
 """
mutable struct ManyBodySolver
    mesh     ::Mesh
    beta     ::Float64
    U        ::Float64
    UX       ::Float64
    n        ::Float64
    sfc_tol  ::Float64
    maxiter  ::Int64
    U_maxiter::Int64
    mix      ::Float64
    verbose  ::Bool
    mu       ::Float64
    gkio     ::Array{ComplexF32,4}
    grit     ::Array{ComplexF32,4}
    ckio     ::Array{ComplexF32,4}
    V        ::Array{ComplexF32,4}
    sigma    ::Array{ComplexF32,4}
    iw       ::Array{ComplexF32,4}
    iv       ::Array{ComplexF32,4}
    ek       ::Array{Float32,4}
    dim      ::Int64
    nk       ::Int64
end

"""Initiarize function"""
function make_ManyBodySolver(
        mesh      ::Mesh,
        beta      ::Float64,
        ek        ::Array{Float32, 4},
        U         ::Float64,
        mu         ::Float64,
        sigma_init::Array{ComplexF32,4};
        sfc_tol   ::Float64=1e-4,
        maxiter   ::Int64  =100,
        U_maxiter ::Int64  =10,
        mix       ::Float64=0.2,
        verbose   ::Bool   =true
        )::ManyBodySolver
    
        n::Float64 = 0.0
    
        gkio  = Array{ComplexF32}(undef, mesh.fnw,   nk1, nk2, nk3)
        grit  = Array{ComplexF32}(undef, mesh.fntau, nk1, nk2, nk3)
        ckio  = Array{ComplexF32}(undef, mesh.bnw,   nk1, nk2, nk3)
        V     = Array{ComplexF32}(undef, mesh.bntau, nk1, nk2, nk3)
        sigma = sigma_init
    
        iw, iv = get_iw_iv(mesh)
        iw = reshape(iw, mesh.fnw, 1, 1, 1)
        iv = reshape(iv, mesh.bnw, 1, 1, 1)
        solver = ManyBodySolver(mesh, beta, U, 0.0, n, sfc_tol, maxiter, U_maxiter, mix, verbose, mu, gkio, grit, ckio, V, sigma, iw, iv, ek, dim, nk)
        solver.n = calc_electron_density(solver, mu)
    
        solver.mu = mu_calc(solver)
        solver.gkio .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.sigma)
        solver.grit .= kw_to_rtau(solver.gkio, 'F', solver.mesh)
        solver.ckio .= rtau_to_kw(solver.grit .* reverse(solver.grit, dims=1), 'B', solver.mesh)
        solver.UX = solver.U * maximum(abs, solver.ckio)

        return solver
end

#%%%%%%%%%%% Loop solving instance
function solve!(solver::ManyBodySolver, comm)
    """ ManyBodySolver.solve() executes FLEX loop until convergence """
    # check whether U < U_crit! Otherwise, U needs to be renormalized.
    old_U = solver.U
    #println("U, X, = $(solver.U), $(maximum(abs, solver.ckio))")
    if solver.UX >= 1
        println("U * max(X) = $(solver.UX). U Renormalization Starting")
        U_renormalization(solver, comm)
        println("New U = $(solver.U)")
        if (solver.U / old_U < 0.9)
            println("-----------------------------------------------")
            println("U is too large! Paramagnetic Phase Unavoidable!")
            println("-----------------------------------------------")
            exit()
        end
    end
            
    println("Beginning Self-Consistent Loop")
    if solver.verbose println("---------------\nIter  | Error\n---------------") end
    # perform loop until convergence is reached:
    sigma_old = copy(solver.sigma)
    #println("U, X, = $(solver.U), $(maximum(abs, solver.ckio))")
    for it in 1:solver.maxiter
        #println("U, X, = $(solver.U), $(maximum(abs, solver.ckio))")
        loop!(solver)
        if solver.UX >= 1 
            println("Divergence in interaction found. Paramagnetic phase entered. Code exiting")
            exit()
        end
        
        # check whether solution is converged.
        sfc_check = sum(abs.(solver.sigma-sigma_old))/sum(abs.(solver.sigma))

        if solver.verbose
            println(it, '\t', sfc_check)
        end
        if sfc_check < solver.sfc_tol
            println("FLEX loop converged at desired accuracy")
            break
        end
        sigma_old .= copy(solver.sigma)
    end
    solver.grit .= kw_to_rtau(solver.gkio, 'F', solver.mesh)
    solver.ckio .= rtau_to_kw(solver.grit .* reverse(solver.grit, dims=1), 'B', solver.mesh)
end
    
function FLEX_loop!(solver::ManyBodySolver, comm)
    """ FLEX loop """
    gkio_old = copy(solver.gkio)
    
    V_calc(solver)
    @inbounds @simd for i in eachindex(solver.V)
        solver.grit[i] = solver.V[i] * solver.grit[i]
    end
    rtau_to_kw!(solver.sigma, solver.grit, 'F', solver.mesh)
        
    #local_sum = sum(real, solver.gkio)
    #global_sum = MPI.Allreduce(local_sum, +, comm)
    #total_size = solver.nk * solver.mesh.fnw
    #solver.mu = Float32(global_sum / total_size)
    solver.mu = mu_calc(solver)
    solver.gkio .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.sigma)
    
    @inbounds @simd for i in eachindex(solver.gkio)
        solver.gkio[i] = solver.mix*solver.gkio[i] + (1-solver.mix)*gkio_old[i]
    end
        
    kw_to_rtau!(solver.grit, solver.gkio, 'F', solver.mesh)        
    rtau_to_kw!(solver.ckio, solver.grit .* reverse(solver.grit, dims=1), 'B', solver.mesh)        
end

function DMFT_Sigma(solver::ManyBodySolver)
    proj = ones(ComplexF32, nx, ny, nz)
    g_loc = project_kernel(proj, solver.gkio)
    gt = wn_to_tau(solver.mesh, Fermionic(), g_loc)
    sigma_t = solver.U^2 * gt.^3
    sigma_w = tau_to_wn(solver.mesh, Fermionic(), sigma_t)
    return reshape(sigma_w, :, 1, 1, 1)
end

function FLEX_local_Sigma(solver::ManyBodySolver)
    proj = ones(ComplexF32, nx, ny, nz)
    g_loc = project_kernel(proj, solver.gkio)
    gt = wn_to_tau(solver.mesh, Fermionic(), g_loc)
    xt = gt .* reverse(gt)
    xw = tau_to_wn(solver.mesh, Bosonic(), xt)
    vw = similar(xw)
    V_FLEX!(solver.U, xw, vw)
    vt = wn_to_tau(solver.mesh, Bosonic(), vw)
    sigmat = vt .* gt
    sigmaw = tau_to_wn(solver.mesh, Fermionic(), sigmat)
    return sigmaw
end

function FLEX_nonlocal_Sigma(solver::ManyBodySolver, local_sigma)
    return solver.sigma .- reshape(local_sigma, :, 1, 1, 1)
end

function FLEX_sigma!(solver::ManyBodySolver)
    # G Creation
    kw_to_rtau!(solver.grit, solver.gkio, 'F', solver.mesh)
    rtau_to_kw!(solver.ckio, solver.grit .* reverse(solver.grit, dims=1), 'B', solver.mesh)
    solver.UX = solver.U * maximum(abs, solver.ckio)
    V_calc(solver)
    rtau_to_kw!(solver.sigma, solver.V .* solver.gkio, 'F', solver.mesh)
end


function combine_sigmas!(solver::ManyBodySolver, sigma_local)
    solver.sigma .= solver.sigma .+ reshape(sigma_local, :, 1, 1, 1)
end


function DFMT_FLEX!(solver::ManyBodySolver)
    # Loop 1
    #println("1) U, X, = $(solver.U), $(maximum(abs, solver.ckio))")
    sigma_w = DMFT_Sigma(solver)
    combine_sigmas!(solver, sigma_w)
    solver.mu = mu_calc(solver)
    solver.gkio .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.sigma)
    #println("2) U, X, = $(solver.U), $(maximum(abs, solver.ckio))")

    FLEX_sigma!(solver)

    # Loop 2
    sigma_w .= FLEX_local_Sigma(solver)
    solver.sigma .= FLEX_nonlocal_Sigma(solver, sigma_w)
    combine_sigmas!(solver, sigma_w)
    solver.mu = mu_calc(solver)
    solver.gkio .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.sigma)
end

function loop!(solver::ManyBodySolver)
    if interaction == "FLEX"
        FLEX_loop!(solver, 1)
    elseif interaction == "DMFT+FLEX" || interaction == "FLEX+DMFT"
        DFMT_FLEX!(solver)
    else
        println("Interaction $interaction not available")
        exit()
    end
end


#%%%%%%%%%%% U renormalization loop instance
function U_renormalization(solver::ManyBodySolver, comm)
    """ Loop for renormalizing U if Stoner enhancement U*max{chi0} >= 1. """
    println("WARNING: U is too large and the spin susceptibility denominator will diverge/turn unphysical!")
    println("Initiate U renormalization loop.")
    
    # save old U for later
    U_old::Float64 = solver.U
    # renormalization loop may run infinitely! Insert break condition after U_it_max steps
    U_it::Int64 = 0

    UX = solver.UX 
    prev_U = 0.0
    U_diff = 1.0
    
    while UX >= 1.0
        U_it += 1
        # reset U
        solver.U = U_old / (solver.UX + 0.01)
        #solver.U = U_old
        
        # remormalize U such that U*chi0 < 1
        println(U_it, '\t', solver.U, '\t', U_old)
        
        # perform one shot FLEX loop
        loop!(solver)
        solver.UX *= U_old / solver.U
        
        
        # break condition for too many steps
        if U_it == solver.U_maxiter || U_diff < 1e-4
            println("U_diff = $U_diff")
            println("Iterations = $U_it")
            println("Maximum Iterations = $(solver.U_maxiter)")
            break
        end
        U_diff = abs(solver.U - prev_U)
        prev_U = solver.U
    end
    solver.UX = solver.U * maximum(abs, solver.ckio)
    println("Leaving U renormalization...")
end

function V_FLEX!(U, ckio, out)
    U2 = U^2
    @inbounds @simd for i in eachindex(ckio)
        x = ckio[i]
        term1 = (1.5*U^2) * x / (1 - U*x)
        term2 = (0.5*U^2) * x / (1 + U*x)
        out[i] = term1 + term2 - U2 * x
    end
end

function V_calc(solver::ManyBodySolver)
    #println("U, X, = $(solver.U), $(maximum(abs, solver.ckio))")
    maxval = maximum(abs.(solver.ckio))*solver.U
    if maxval >= 1
        error("U*max(chi0) = $(maxval) >= 1! Paramagnetic phase is left and calculations will turn unstable!")
    end

    V_FLEX!(solver.U, solver.ckio, solver.ckio)
    # Constant Hartree Term V ~ U needs to be treated extra, since they cannot be modeled by the IR basis.
    # In the single-band case, the Hartree term can be absorbed into the chemical potential.

    #solver.V .= kw_to_rtau(solver.ckio, 'B', solver.mesh)
    kw_to_rtau!(solver.V, solver.ckio, 'B', solver.mesh)
end

#%%%%%%%%%%% Setting chemical potential mu
function calc_electron_density(solver::ManyBodySolver,mu::Float64)::Float64
    """ Calculate electron density from Green function """
    solver.gkio .= 1.0 ./ (solver.iw .- (solver.ek .- mu) .- solver.sigma)
    gio = dropdims(sum(solver.gkio,dims=(2,3)),dims=(2,3))/solver.mesh.nk

    g_l = fit(solver.mesh.IR_basis_set.smpl_wn_f,gio, dim=1)
    g_tau0 = dot(solver.mesh.IR_basis_set.basis_f.u(0), g_l)

    n  = 1.0 + real(g_tau0)
    n  = 2.0 * n #for spin
end

function mu_calc(solver::ManyBodySolver)::Float64
    n = Float64(solver.n)
    """ Find chemical potential for a given filling n0 via brent's root finding algorithm """
    f  = x -> calc_electron_density(solver,Float64(x)) - n

    mu = find_zero(f, (3*Float64(minimum(solver.ek)), 3*Float64(maximum(solver.ek))), Roots.Brent()) 
    return mu
end

function get_local_k_slice(rank, nprocs, total_kz)
    base = div(total_kz, nprocs)
    rem = mod(total_kz, nprocs)

    start_kz = rank * base + min(rank, rem) + 1
    end_kz = start_kz + base - 1
    if rank < rem
        end_kz += 1
    end
    return start_kz:end_kz
end
# initialize calculation

function main()

    mpi_test()

    comm = 1

    band = Bands()
    ek = fill_energy_mesh(band)
    mesh = IR_Mesh(ek)

    sigma_init = zeros(ComplexF32,(mesh.fnw, nk1, nk2, nk3))
    verbose = cfg.verbosity == "high" 
    solver = make_ManyBodySolver(mesh, beta, ek, U, mu, sigma_init, sfc_tol=sfc_tol, maxiter=maxiter, U_maxiter=U_maxiter, mix=mix, verbose=verbose)

    # perform FLEX loop
    if scf
        solve!(solver, comm)
    end
    println("New mu=$(solver.mu)")

    solver.grit = Array{ComplexF32,4}(undef, 0, 0, 0, 0)
    solver.V = Array{ComplexF32,4}(undef, 0, 0, 0, 0)
    V = similar(solver.ckio)

    V_FLEX!(solver.U, solver.ckio, V)

    println("Max Self-Energy: $(maximum(abs.(solver.sigma)))")
    println("Min Self-Energy: $(minimum(abs.(solver.sigma)))")
    println("Max Vertex: $(maximum(abs.(V)))")
    println("Min Vertex: $(minimum(abs.(V)))")
    println("Max Chi: $(maximum(abs.(solver.ckio)))")
    println("Min Chi: $(minimum(abs.(solver.ckio)))")
    #if nz > 1
    #    sitp = interpolate((1:mesh.fnw, 1:nx, 1:ny, 1:nz), Sigma, Gridded(Linear()))
    #    vitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:nz), V, Gridded(Linear()))
    #    xitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:nz), X, Gridded(Linear()))
    #else
    #    Sigma = repeat(Sigma, 1, 1, 1, 2)
    #    V = repeat(V, 1, 1, 1, 2)
    #    X = repeat(X, 1, 1, 1, 2)
    #    sitp = interpolate((1:mesh.fnw, 1:nx, 1:ny, 1:2), Sigma, Gridded(Linear()))
    #    vitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:2), V, Gridded(Linear()))
    #    xitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:2), X, Gridded(Linear()))
    #end

    BZ_in = BZ
    kmesh = cfg.k_mesh
    if dim == 2
        kmesh = kmesh[1:end-1]
        BZ_in = BZ_in[1:end-1, 1:end-1]
    end

    # Centers points correctly, so they go from (-pi,pi) to (pi,pi) instead of the current (0,0) to (2pi,2pi). Important for saving
    for i in 1:mesh.fnw
        solver.sigma[i, :, :, :] .= fftshift(solver.sigma[i, :, :, :])
    end
    for i in 1:mesh.bnw
        V[i, :, :, :] .= fftshift(V[i, :, :, :])
        solver.ckio[i, :, :, :] .= fftshift(solver.ckio[i, :, :, :])
    end

    save_field!(outdir * prefix * "_self_energy." * filetype, solver.sigma, BZ_in, kmesh, imag.(solver.iw))
    save_field!(outdir * prefix * "_vertex." * filetype, V, BZ_in, kmesh, imag.(solver.iv))
    save_field!(outdir * prefix * "_chi." * filetype, solver.ckio, BZ_in, kmesh, imag.(solver.iv))

end

function mpi_test()
    MPI.Init()
    comm = MPI.COMM_WORLD

    #rank = MPI.Comm_rank(comm)
    #nprocs = MPI.Comm_size(comm)
    #println("comm = $comm rank = $rank nprocs = $nprocs")

    band = Bands()
    ek = fill_energy_mesh(band)
    @show size(ek)
    mesh = IR_Mesh(ek)
    iw, iv = get_iw_iv(mesh)
    @show size(iw)

    dims = (nk1, nk2, nk3)
    pen = Pencil(dims, comm)
    transform = Transforms.FFT()
    farg = Val(false)
    plan_fnw = PencilFFTPlan(pen, transform, extra_dims = (mesh.fnw,), permute_dims = farg)
    plan_bnw = PencilFFTPlan(pen, transform, extra_dims = (mesh.bnw,), permute_dims = farg)
    plan_fntau = PencilFFTPlan(pen, transform, extra_dims = (mesh.fntau,), permute_dims = farg)
    plan_bntau = PencilFFTPlan(pen, transform, extra_dims = (mesh.bntau,), permute_dims = farg)
    @show size(plan_fnw)
    @show size(plan_bnw)
    @show size(plan_fntau)
    @show size(plan_bntau)

    smpl_tau_F = mesh.IR_basis_set.smpl_tau_f
    smpl_wn_F = mesh.IR_basis_set.smpl_wn_f
    smpl_tau_B = mesh.IR_basis_set.smpl_tau_b
    smpl_wn_B = mesh.IR_basis_set.smpl_wn_b
    G = allocate_input(plan_fnw)
    G_rw = allocate_output(plan_fnw)
    #temp_frw = allocate_output(plan_fnw)
    temp = allocate_output(plan_fntau)
    temp_frt = allocate_input(plan_fntau)
    temp_brw = allocate_input(plan_bnw)
    G_rt = allocate_output(plan_fntau)
    G1_rt = allocate_input(plan_fntau)
    X_rt = allocate_input(plan_bntau)
    X = allocate_output(plan_bnw)
    @show size(G)
    @show size(temp_frt)
    @show size(temp_brw)
    @show size(G_rt)
    @show size(X_rt)
    @show size(X)

    println("0")
    G .= 1 ./ (reshape(iw, 1, 1, 1, :) .- dropdims(ek; dims=1) .+ mu)
    println("1")
    mul!(G_rw, plan_fnw, G)
    temp_frw = similar(G)
    #G_rw = plan_fnw * G
    println("2")
    ldiv!(temp_frw, plan_fnw, G_rw)
    println("3")
    u = plan_fnw / G_rw
    println("4")
    println("done")
    println("G->l")
    obj_l_F = fit(smpl_wn_F, parent(G), dim=4)
    @show size(obj_l_F)
    println("l->t")
    evaluate!(temp_frt, smpl_tau_F, obj_l_F, dim=4)
    println("kt->rt")
    mul!(G_rt, plan_fntau, temp_frt)
    #temp = similar(temp_frt)
    #G1_rt .= G_rt
    #mul!(temp, adjoint(plan_fntau), G1_rt)  # now w ≈ u
    #ldiv!(temp, plan_fntau, G1_rt)  # now w ≈ u
    println("rt->rt")
    X_rt = G_rt .* reverse(G_rt, dims=(4))
    println("rt->l")
    obj_l_B = fit(smpl_tau_B, parent(X_rt), dim=4)
    @show size(obj_l_B)
    println("l->rw")
    evaluate!(temp_brw, smpl_wn_B, obj_l_B, dim=4)
    #println("rw->kw")
    #temp = plan_bnw / temp_brw
    println("rw->kw")
    ldiv!(X, plan_bnw, temp_brw)  # now w ≈ u
    @show size(X)

    println("Max X: $(maximum(abs, X))")

    println("Rank $rank: completed test_flex")

    MPI.Barrier(comm)
    MPI.Finalize()
    return
end

end # module
