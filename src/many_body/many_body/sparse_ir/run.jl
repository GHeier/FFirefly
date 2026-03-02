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
using DelimitedFiles
#using Interpolations
include("../../../objects/mesh.jl")
using .IRMesh
using Firefly
cfg = Firefly.Config

outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype

scf = cfg.self_consistent 
scf_tol = 1e-4
mix = cfg.mixing
mu_from_n = cfg.mu_from_n

nx, ny, nz = cfg.k_mesh
nqx, nqy, nqz = cfg.q_mesh
dim = cfg.dimension
if dim == 2
    nz = 1
    nqz = 1
end
nk = nx * ny * nz
nbnd = cfg.nbnd
nbnd = cfg.nbnd
mu = cfg.fermi_energy
n = cfg.num_electrons
U = cfg.U0
BZ = cfg.brillouin_zone
beta = 1 / cfg.Temperature


### System parameters
#n    = 0.85   # electron filling, here per spin per lattice site (n=1: half filling)

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
        ek[i, j, k] = band(kvec)
    end
    return ek
end

function fill_energy_mesh(band)
    ek = Array{Float32}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i - 1, j - 1, k - 1, nx, ny, nz)
        ek[i, j, k] = band(kvec)
    end
    return ek
end

function get_energy_min_max(ek)
    return minimum(ek), maximum(ek)
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
    Gkw      ::Array{ComplexF32, 4}
    Grt      ::Array{ComplexF32, 4}
    Xkw      ::Array{ComplexF32, 4}
    Vrt      ::Array{ComplexF32, 4}
    Ekw      ::Array{ComplexF32, 4}
    iw       ::Array{ComplexF32,4}
    iv       ::Array{ComplexF32,4}
    ek       ::Array{Float32, 4}
    dim      ::Int64
    nk       ::Int64
end


"""Initiarize function"""
function make_ManyBodySolver(
        mesh      ::Mesh,
        beta      ::Float64,
        ek        ::Array{Float32, 3},
        U         ::Float64,
        mu        ::Float64,
        n         ::Float64,
        sigma_init::Array{ComplexF32, 4};
        sfc_tol   ::Float64=1e-4,
        maxiter   ::Int64  =100,
        U_maxiter ::Int64  =10,
        mix       ::Float64=0.2,
        verbose   ::Bool   =true
        )::ManyBodySolver
    
        #n::Float64 = 0.0
    
        #Gkw = Array{Array{ComplexF32,4},4}(undef, mesh.fnw, 
        Gkw  = zeros(ComplexF32, mesh.fnw, nk1, nk2, nk3)
        Grt  = zeros(ComplexF32, mesh.fntau, nk1, nk2, nk3)
        Xkw  = zeros(ComplexF32, mesh.bnw, nk1, nk2, nk3)
        V    = zeros(ComplexF32, mesh.bntau, nk1, nk2, nk3)
        Ekw = sigma_init
    
        iw, iv = get_iw_iv(mesh)
        iw = reshape(iw, mesh.fnw, 1, 1, 1)
        iv = reshape(iv, mesh.bnw, 1, 1, 1)
        ek_new = reshape(ek, 1, nx, ny, nz)
        solver = ManyBodySolver(mesh, beta, U, 0.0, n, sfc_tol, maxiter, U_maxiter, mix, verbose, mu, Gkw, Grt, Xkw, V, Ekw, iw, iv, ek_new, dim, nk)
        if mu_from_n
            solver.n = n
            println("Initial mu = $mu")
            solver.mu = mu_calc(solver)
            println("New mu = $(solver.mu)")
        else
            solver.mu = mu
            println("Initial n = $n")
            solver.n = calc_electron_density(solver, mu)
            println("New n = $(solver.n)")
        end
    
        solver.Gkw .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.Ekw)
        solver.Grt .= kw_to_rtau(solver.Gkw, 'F', solver.mesh)
        solver.Xkw .= rtau_to_kw(solver.Grt .* reverse(solver.Grt, dims=1), 'B', solver.mesh)
        solver.UX = solver.U * maximum(abs, solver.Xkw)
        println("Initial Max Chi = $(maximum(abs, solver.Xkw))")
        println("Initial UX = $(solver.UX)")

        return solver
end

#%%%%%%%%%%% Loop solving instance
function solve!(S::ManyBodySolver, comm)
    """ ManyBodyS.solve() executes FLEX loop until convergence """
    # check whether U < U_crit! Otherwise, U needs to be renormalized.
    old_U = S.U
    #println("U, X, = $(S.U), $(maximum(abs, S.ckio))")
    if S.UX >= 1 && occursin("FLEX", interaction)
        if scf
            println("U * max(X) = $(S.UX). U Renormalization Starting")
            U_renormalization(S, comm)
            println("New U = $(S.U)")
        else
            println("U * max(X) = $(S.UX). Calculations Impossible")
            exit(1)
        end
        if (S.U / old_U < 0.9)
            println("-----------------------------------------------")
            println("U is too large! Paramagnetic Phase Unavoidable!")
            println("Continuing with reduced U...")
            println("-----------------------------------------------")
        end
        # Reset U to original value after renormalization (matching fulltest.jl behavior)
        S.U = old_U
    end
            
    println("Beginning Self-Consistent Loop")
    if S.verbose println("---------------\nIter  | Error\n---------------") end
    # perform loop until convergence is reached:
    sigma_old = copy(S.Ekw)
    #println("U, X, = $(S.U), $(maximum(abs, S.ckio))")
    for it in 1:S.maxiter
        loop!(S)

        # check whether solution is converged.
        sfc_check = sum(abs.(S.Ekw-sigma_old))/sum(abs.(S.Ekw))

        if S.verbose
            println(it, '\t', sfc_check)
        end
        if sfc_check < S.sfc_tol
            println("FLEX loop converged at desired accuracy")
            break
        end
        sigma_old .= copy(S.Ekw)
    end
    #S.Grt .= kw_to_rtau(S.Gkw, 'F', S.mesh)
    #S.Xkw .= rtau_to_kw(S.Grt .* reverse(S.Grt, dims=1), 'B', S.mesh)
end
    
function FLEX_loop!(S::ManyBodySolver, comm)
    """ FLEX loop """
    gkio_old = copy(S.Gkw)
    
    V_calc(S)
    S.Grt .= S.Vrt .* S.Grt
    rtau_to_kw!(S.Ekw, S.Grt, 'F', S.mesh)
    #@inbounds @simd for i in eachindex(S.V)
    #    S.grit[i] = S.V[i] * S.grit[i]
    #end
        
    #local_sum = sum(real, S.gkio)
    #global_sum = MPI.Allreduce(local_sum, +, comm)
    #total_size = S.nk * S.mesh.fnw
    #S.mu = Float32(global_sum / total_size)
    S.mu = mu_calc(S)
    S.Gkw .= 1.0 ./ (S.iw .- (S.ek .- S.mu) .- S.Ekw)
    
    S.Gkw .= S.mix.*S.Gkw .+ (1-S.mix).*gkio_old
        
    kw_to_rtau!(S.Grt, S.Gkw, 'F', S.mesh)        
    S.Xkw .= rtau_to_kw(S.Grt .* reverse(S.Grt, dims=1), 'B', S.mesh)
end

function DMFT_Sigma(solver::ManyBodySolver)
    proj = ones(ComplexF32, nx, ny, nz)
    g_loc = project_kernel(proj, solver.Gkw)
    gt = wn_to_tau(solver.mesh, Fermionic(), g_loc)
    sigma_t = solver.U^2 * gt.^3
    sigma_w = tau_to_wn(solver.mesh, Fermionic(), sigma_t)
    return reshape(sigma_w, :, 1, 1, 1)
end

function FLEX_local_Sigma(solver::ManyBodySolver)
    proj = ones(ComplexF32, nx, ny, nz)
    xw = zeros(ComplexF32, solver.mesh.bnw)
    gt = zeros(ComplexF32, solver.mesh.fnw)
    g_loc = project_kernel(proj, solver.Gkw)
    gt .= wn_to_tau(solver.mesh, Fermionic(), g_loc)
    xt = gt .* reverse(gt)
    xw .+= tau_to_wn(solver.mesh, Bosonic(), xt)
    vw = similar(xw)
    V_FLEX!(solver.U, xw, vw)
    vt = wn_to_tau(solver.mesh, Bosonic(), vw)
    sigmaw = zeros(ComplexF32, nbnd, solver.mesh.fnw)
    sigmat = vt .* gt
    sigmaw .= tau_to_wn(solver.mesh, Fermionic(), sigmat)
    return sigmaw
end

function FLEX_nonlocal_Sigma(solver::ManyBodySolver, local_sigma)
    return solver.Ekw .- reshape(local_sigma, :, 1, 1, 1)
end

function FLEX_sigma!(solver::ManyBodySolver)
    # G Creation
    kw_to_rtau!(solver.Grt, solver.Gkw, 'F', solver.mesh)
    solver.Xkw .+= rtau_to_kw(solver.Grt .* reverse(solver.Grt, dims=1), 'B', solver.mesh)
    #rtau_to_kw!(solver.Xkw[1], solver.Grt .* reverse(solver.Grt, dims=1), 'B', solver.mesh)
    solver.UX = solver.U * maximum(abs, solver.Xkw)
    V_calc(solver)
    rtau_to_kw!(solver.Ekw, solver.Vrt .* solver.Gkw, 'F', solver.mesh)
end


function combine_sigmas!(solver::ManyBodySolver, sigma_local)
    solver.Ekw .= solver.Ekw .+ reshape(sigma_local, :, 1, 1, 1)
end


function DFMT_FLEX!(solver::ManyBodySolver)
    # Loop 1
    sigma_w = DMFT_Sigma(solver)
    combine_sigmas!(solver, sigma_w)
    solver.mu = mu_calc(solver)
    solver.Gkw .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.Ekw)

    FLEX_sigma!(solver)

    # Loop 2
    sigma_w .= FLEX_local_Sigma(solver)
    solver.Ekw .= FLEX_nonlocal_Sigma(solver, sigma_w)
    combine_sigmas!(solver, sigma_w)
    solver.mu = mu_calc(solver)
    solver.Gkw .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.Ekw)
end

function DMFT_loop!(solver::ManyBodySolver)
    sigma_w = DMFT_Sigma(solver)
    combine_sigmas!(solver, sigma_w)
    solver.mu = mu_calc(solver)
    solver.Gkw .= 1.0 ./ (solver.iw .- (solver.ek .- solver.mu) .- solver.Ekw)
end

function loop!(solver::ManyBodySolver)
    if interaction == "FLEX"
        FLEX_loop!(solver, 1)
    elseif interaction == "DMFT"
        DMFT_loop!(solver)
    elseif interaction == "DMFT+FLEX" || interaction == "FLEX+DMFT"
        DFMT_FLEX!(solver)
    else
        println("Interaction $interaction not available")
        exit()
    end
end


#%%%%%%%%%%% U renormalization loop instance
function U_renormalization(solver::ManyBodySolver, comm)
    """ Loop for renormalizing U if Stoner enhancement U*max{chi0} >= 1.
        Matches fulltest.jl behavior: temporarily reduce U, run loop, reset U. """
    println("WARNING: U is too large and the spin susceptibility denominator will diverge/turn unphysical!")
    println("Initiate U renormalization loop.")

    # save old U for later
    U_old::Float64 = solver.U
    U_it::Int64 = 0

    while U_old * maximum(abs, solver.Xkw) >= 1.0
        U_it += 1

        # Temporarily reduce U so U*chi0 < 1
        solver.U = solver.U / (maximum(abs, solver.Xkw) * solver.U + 0.01)
        println(U_it, '\t', solver.U, '\t', U_old)

        # Perform one FLEX loop iteration
        loop!(solver)

        # Reset U to original value (key difference from before)
        solver.U = U_old

        # Break condition for too many steps
        if U_it == solver.U_maxiter
            println("U renormalization reached breaking point")
            break
        end
    end
    println("Leaving U renormalization...")
end

function V_FLEX!(U, Xkw, out)
    # Triplet channel vertex for self-energy: 1.5*U²*χ_spin + 0.5*U²*χ_charge - U²*χ₀
    out .= (1.5*U^2) .* Xkw ./ (1 .- U.*Xkw) .+ (0.5*U^2) .* Xkw ./ (1 .+ U.*Xkw) - U^2 .* Xkw
end

function V_singlet!(U, Xkw, out)
    # Singlet channel vertex for pairing: 1.5*U²*χ_spin - 0.5*U²*χ_charge
    out .= (1.5*U^2) .* Xkw ./ (1 .- U.*Xkw) .- (0.5*U^2) .* Xkw ./ (1 .+ U.*Xkw)
end

function V_calc(solver::ManyBodySolver)
    #println("U, X, = $(solver.U), $(maximum(abs, solver.ckio))")
    maxval = maximum(abs.(solver.Xkw))*solver.U
    if maxval >= 1
        error("U*max(chi0) = $(maxval) >= 1! Paramagnetic phase is left and calculations will turn unstable!")
    end

    V_FLEX!(solver.U, solver.Xkw, solver.Xkw)
    # Constant Hartree Term V ~ U needs to be treated extra, since they cannot be modeled by the IR basis.
    # In the single-band case, the Hartree term can be absorbed into the chemical potential.

    #solver.V .= kw_to_rtau(solver.ckio, 'B', solver.mesh)
    kw_to_rtau!(solver.Vrt, solver.Xkw, 'B', solver.mesh)
end

#%%%%%%%%%%% Setting chemical potential mu
function calc_electron_density(S::ManyBodySolver,mu::Float64)::Float64
    """ Calculate electron density from Green function """
    S.Gkw .= 1.0 ./ (S.iw .- (S.ek .- mu) .- S.Ekw)
    gio = dropdims(sum(S.Gkw,dims=(2,3,4)),dims=(2,3,4))/S.mesh.nk

    g_l = fit(S.mesh.IR_basis_set.smpl_wn_f,gio, dim=1)
    g_tau0 = dot(S.mesh.IR_basis_set.basis_f.u(0), g_l)

    n_tmp  = 1.0 + real(g_tau0)
    n_tmp  = 2.0 * n_tmp #for spin
    return n_tmp
end

function mu_calc(solver::ManyBodySolver)::Float64
    n = Float64(solver.n)
    """ Find chemical potential for a given filling n0 via brent's root finding algorithm """
    f  = x -> calc_electron_density(solver,Float64(x)) - n


    minval, maxval = get_energy_min_max(solver.ek)
    mu = find_zero(f, (3*Float64(minval), 3*Float64(maxval)), Roots.Brent()) 
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

function find_zero_freq(iw)
    for i in 1:length(iw)
        if imag(iw[i]) * imag(iw[i+1]) < 0
            return i
        end
    end
end

function main()

    #mpi_test2()
    #return

    comm = 1
    println("Constructing Bands")
    band = Bands()

    println("Filling Energy Mesh")
    ek = fill_energy_mesh(band)
    minval, maxval = get_energy_min_max(ek)
    println("Minimum Energy = $(minval)")
    println("Maximum Energy = $(maxval)")
    D = maxval - minval
    mesh = IR_Mesh(1.2*D)

    iw, iv = get_iw_iv(mesh)
    sigma_init = zeros(ComplexF32, mesh.fnw, nk1, nk2, nk3)

    verbose = cfg.verbosity == "high" 
    solver = make_ManyBodySolver(mesh, beta, ek, U, mu, n, sigma_init, sfc_tol=sfc_tol, maxiter=maxiter, U_maxiter=U_maxiter, mix=mix, verbose=verbose)

    X_copy = copy(solver.Xkw)

    Gkw = 1.0 ./ (reshape(iw, mesh.fnw, 1, 1, 1) .- (reshape(ek, 1, nx, ny, nz) .+ solver.mu))
    ind = Int(mesh.fnw / 2)
    G_w0 = sum(Gkw[ind, :, :, :]) / nk
    println("Non-interacting DOS = $(G_w0.im / pi)")

    # perform FLEX loop
    # Always run solve! - maxiter is already set to 1 when scf=false (line 54)
    solve!(solver, comm)
    if !scf
        solver.Xkw .= X_copy
    end
    println("New mu=$(solver.mu)")

    V = similar(solver.Xkw)
    V_singlet = similar(solver.Xkw)

    V_FLEX!(solver.U, solver.Xkw, V)
    V_singlet!(solver.U, solver.Xkw, V_singlet)

    ind = find_zero_freq(iw)
    slope = -(imag(solver.Ekw[ind+1]) - imag(solver.Ekw[ind])) / (imag(iw[ind+1]) - imag(iw[ind]))
    Z = 1 / (1 + slope)
    println("Renormalization: ", 1 + slope)
    println("Quasiparticle Weight: ", Z)

    println("Max Self-Energy: $(maximum(abs.(solver.Ekw)))")
    println("Min Self-Energy: $(minimum(abs.(solver.Ekw)))")
    println("Max Vertex: $(maximum(abs.(V)))")
    println("Min Vertex: $(minimum(abs.(V)))")
    println("Max Chi: $(maximum(abs.(solver.Xkw)))")
    println("Min Chi: $(minimum(abs.(solver.Xkw)))")

    BZ_in = BZ
    kmesh = cfg.k_mesh
    if dim == 2
        kmesh = kmesh[1:end-1]
        BZ_in = BZ_in[1:end-1, 1:end-1]
    end


    G_w0 = sum(solver.Gkw[ind, :, :, :]) / nk
    println("Interacting DOS = $(G_w0.im / pi)")

    iw, iv = reshape(solver.iw, mesh.fnw), reshape(solver.iv, mesh.bnw)
    save_data!(outdir * prefix * "_self_energy." * filetype, solver.Ekw, kmesh, BZ_in, w_points=imag.(iw), centered=false)
    println("Saving to ", outdir * prefix * "_self_energy.h5")

    save_data!(outdir * prefix * "_vertex." * filetype, V, kmesh, BZ_in, w_points=imag.(iv), centered=false)
    println("Saving to ", outdir * prefix * "_vertex.h5")

    save_data!(outdir * prefix * "_vertex_singlet." * filetype, V_singlet, kmesh, BZ_in, w_points=imag.(iv), centered=false)
    println("Saving to ", outdir * prefix * "_vertex_singlet.h5")

    save_data!(outdir * prefix * "_chi." * filetype, solver.Xkw, kmesh, BZ_in, w_points=imag.(iv), centered=false)
    println("Saving to ", outdir * prefix * "_chi.h5")

    return maximum(abs.(solver.Xkw))
end

function mpi_test()
    MPI.Init()
    comm = MPI.COMM_WORLD

    #rank = MPI.Comm_rank(comm)
    #nprocs = MPI.Comm_size(comm)
    #println("comm = $comm rank = $rank nprocs = $nprocs")

    band = Bands()
    ek = fill_energy_mesh(band)
    mesh = IR_Mesh(ek)
    iw, iv = get_iw_iv(mesh)

    dims = (nk1, nk2, nk3)
    pen = Pencil(dims, comm)
    transform = Transforms.FFT()
    farg = Val(false)
    plan_fnw = PencilFFTPlan(pen, transform, extra_dims = (mesh.fnw,), permute_dims = farg)
    plan_bnw = PencilFFTPlan(pen, transform, extra_dims = (mesh.bnw,), permute_dims = farg)
    plan_fntau = PencilFFTPlan(pen, transform, extra_dims = (mesh.fntau,), permute_dims = farg)

    smpl_tau_F = mesh.IR_basis_set.smpl_tau_f
    smpl_wn_F = mesh.IR_basis_set.smpl_wn_f
    smpl_tau_B = mesh.IR_basis_set.smpl_tau_b
    smpl_wn_B = mesh.IR_basis_set.smpl_wn_b

    G = allocate_input(plan_fnw)
    temp_frt = allocate_input(plan_fntau)
    X = allocate_output(plan_bnw)
    G_rt = allocate_output(plan_fntau)
    X_kt = allocate_input(plan_fntau)

    G .= 1 ./ (reshape(iw, 1, 1, 1, :) .- dropdims(ek; dims=1) .+ mu)
    obj_l_F = fit(smpl_wn_F, parent(G), dim=4)
    evaluate!(temp_frt, smpl_tau_F, obj_l_F, dim=4)
    mul!(G_rt, plan_fntau, temp_frt)
    G_rt .= G_rt .* reverse(G_rt, dims=(4))
    ldiv!(X_kt, plan_fntau, G_rt)
    X_kt ./= (nk1 * nk2 * nk3)
    obj_l_B = fit(smpl_tau_B, parent(X_kt), dim=4)
    evaluate!(X, smpl_wn_B, obj_l_B, dim=4)
    println("Max X: $(maximum(abs, X))")

    println("Rank $rank: completed test_flex")


    MPI.Barrier(comm)
    MPI.Finalize()
    return
end

function mpi_test2()
    MPI.Init()
    comm = MPI.COMM_WORLD

    #rank = MPI.Comm_rank(comm)
    #nprocs = MPI.Comm_size(comm)
    #println("comm = $comm rank = $rank nprocs = $nprocs")

    band = Bands()
    ek = fill_energy_mesh(band)
    dims = (nk1, nk2, nk3)
    pen = Pencil(dims, comm)
    mesh = IR_Mesh(ek, pen)
    iw, iv = get_iw_iv(mesh)

    println("allo test")
    G = allocate_input(mesh.plan_fnw)
    temp_frt = allocate_input(mesh.plan_tau)
    X = allocate_output(mesh.plan_bnw)
    G_rt = allocate_output(mesh.plan_tau)
    X_kt = allocate_input(mesh.plan_tau)
    println("allo test done")

    println("1")
    G .= 1 ./ (reshape(iw, 1, 1, 1, :) .- dropdims(ek; dims=1) .+ mu)
    println("2")
    obj_l_F = fit(mesh.smpl_wn_F, parent(G), dim=4)
    println("3")
    evaluate!(temp_frt, mesh.smpl_tau_F, mesh.obj_l_F, dim=4)
    println("4")
    mul!(G_rt, mesh.plan_tau, temp_frt)
    println("5")
    G_rt .= G_rt .* reverse(G_rt, dims=(4))
    ldiv!(X_kt, plan_tau, G_rt)
    X_kt ./= (nk1 * nk2 * nk3)
    obj_l_B = fit(mesh.smpl_tau_B, parent(X_kt), dim=4)
    evaluate!(X, mesh.smpl_wn_B, mesh.obj_l_B, dim=4)
    println("Max X: $(maximum(abs, X))")

    println("Rank $rank: completed test_flex")


    MPI.Barrier(comm)
    MPI.Finalize()
    return
end



function run()
    # Main function call goes here
    max_chi = main()
    return max_chi
end

if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



