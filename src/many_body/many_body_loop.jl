module ManyBodyLoop

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

if cfg.interaction != "FLEX"
    println("Code written for FLEX interactions")
    exit()
end

#function main()
#    band = Bands()
#    ek = fill_energy_mesh(band)
#    mesh = IR_Mesh(ek)
#    iw, iv = get_iw_iv(mesh)
#    fnw, bnw, fntau, bntau = mesh.fnw, mesh.bnw, mesh.fntau, mesh.bntau
#    iw = reshape(iw, fnw, 1, 1, 1)
#    # Preallocate all arrays
#    Sigma = zeros(ComplexF32, fnw, nx, ny, nz)
#    G     = zeros(ComplexF32, fnw, nx, ny, nz)
#    G_rt  = Array{ComplexF32}(undef, fntau, nx, ny, nz)
#    X     = Array{ComplexF32}(undef, bnw, nx, ny, nz)
#    X_rt  = Array{ComplexF32}(undef, bntau, nx, ny, nz)
#    V     = Array{ComplexF32}(undef, bnw, nx, ny, nz)
#
#    solve!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
#
#    println("Max Self-Energy: $(maximum(abs.(Sigma)))")
#    println("Min Self-Energy: $(minimum(abs.(Sigma)))")
#    println("Max Vertex: $(maximum(abs.(V)))")
#    println("Min Vertex: $(minimum(abs.(V)))")
#    println("Max Chi: $(maximum(abs.(X)))")
#    println("Min Chi: $(minimum(abs.(X)))")
#
#    #if nz > 1
#    #    sitp = interpolate((1:mesh.fnw, 1:nx, 1:ny, 1:nz), Sigma, Gridded(Linear()))
#    #    vitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:nz), V, Gridded(Linear()))
#    #    xitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:nz), X, Gridded(Linear()))
#    #else
#    #    Sigma = repeat(Sigma, 1, 1, 1, 2)
#    #    V = repeat(V, 1, 1, 1, 2)
#    #    X = repeat(X, 1, 1, 1, 2)
#    #    sitp = interpolate((1:mesh.fnw, 1:nx, 1:ny, 1:2), Sigma, Gridded(Linear()))
#    #    vitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:2), V, Gridded(Linear()))
#    #    xitp = interpolate((1:mesh.bnw, 1:nx, 1:ny, 1:2), X, Gridded(Linear()))
#    #end
#
#    BZ_in = BZ
#    kmesh = cfg.k_mesh
#    if dim == 2
#        kmesh = kmesh[1:end-1]
#        BZ_in = BZ_in[1:end-1, 1:end-1]
#    end
#
#    # Centers points correctly, so they go from (-pi,pi) to (pi,pi) instead of the current (0,0) to (2pi,2pi). Important for saving
#    for i in 1:mesh.fnw
#        Sigma[i, :, :, :] .= fftshift(Sigma[i, :, :, :])
#    end
#    for i in 1:mesh.bnw
#        V[i, :, :, :] .= fftshift(V[i, :, :, :])
#        X[i, :, :, :] .= fftshift(X[i, :, :, :])
#    end
#
#    save_field!(outdir * prefix * "_self_energy." * filetype, Sigma, BZ_in, kmesh, imag.(iw))
#    save_field!(outdir * prefix * "_vertex." * filetype, V, BZ_in, kmesh, imag.(iv))
#    save_field!(outdir * prefix * "_chi." * filetype, X, BZ_in, kmesh, imag.(iv))
#end
#
#function solve!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
#    max_iters = 30
#    if !scf
#        max_iters = 1
#    end
#
#    G .= 1 ./ (iw .- ek .+ mu .- Sigma)
#    G_rt .= kw_to_rtau(G, 'F', mesh)
#    X_rt .= G_rt .* reverse(G_rt, dims=1)
#    X .= rtau_to_kw(X_rt, 'B', mesh)
#
#    max_val = maximum(abs, X) * U
#    if (max_val >= 1)
#        println("U*X = $(max_val)")
#        if scf
#            U_renormalization(U, mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
#            println("New U = $U")
#        else
#            println("Diverging Interaction, cannot be handled without scf U renormalization. Exiting")
#            exit(1)
#        end
#    end
#
#    for _ in 1:max_iters
#        sigma_old = copy(Sigma)
#        loop!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V, U)
#        scf_check = sum(abs.(Sigma-sigma_old))/sum(abs.(Sigma))
#        println("Max(X) = $(maximum(abs.(X))), err = $(scf_check)")
#        if scf_check < scf_tol
#            println("SCF Condition met")
#            break
#        end
#    end
#end
#
#function loop!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V, U)
#    G_old = copy(G)
#
#    V .= (3/2) .* (U^2 .* X ./ (1 .- U .* X)) .- (1/2) .* (U^2 .* X ./ (1 .+ U .* X))
#    #V .= (U^2 .* X ./ (1 .- U .* X)) .+ (U^3 .* X.^2 ./ (1 .+ U^2 .* X.^2))
#
#    X_rt .= kw_to_rtau(V, 'B', mesh)
#
#    G_rt .= X_rt .* G_rt
#
#    Sigma .= rtau_to_kw(G_rt, 'F', mesh)
#
#    n = calc_electron_density(G, mesh)
#    mu = mu_calc(iw, ek, n, mesh)
#
#    G .= 1 ./ (iw .- ek .+ mu .- Sigma)
#    G .= mix * G .+ (1 - mix) * G_old
#
#    G_rt .= kw_to_rtau(G, 'F', mesh)
#
#    X_rt .= G_rt .* reverse(G_rt, dims=1)
#
#    X .= rtau_to_kw(X_rt, 'B', mesh)
#end
#
#function U_renormalization!(U, mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
#    """ Loop for renormalizing U if Stoner enhancement U*max{chi0} >= 1. """
#    println("WARNING: U is too large and the spin susceptibility denominator will diverge/turn unphysical!")
#    println("Initiate U renormalization loop.")
#
#    # save old U for later
#    U_old::Float64 = U
#    # renormalization loop may run infinitely! Insert break condition after U_it_max steps
#    U_it::Int64 = 0
#    U_maxiter = 50
#
#    println("Iter \t U \t U_old \t max_Chi")
#    while U_old * maximum(abs, X) >= 1.0
#        U_it += 1
#
#        # remormalize U such that U*chi0 < 1
#        U = U / (maximum(abs, X) * U + 0.01)
#
#        # perform one shot FLEX loop
#        loop!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V, U)
#        println(U_it, '\t', round(U, digits = 2), "\t ", U_old, '\t', maximum(abs, X))
#
#        # reset U
#        U = U_old
#
#        # break condition for too many steps
#        if U_it == U_maxiter
#            println("U renormalization reached maximum iterations")
#            break
#        end
#    end
#    println("Leaving U renormalization...")
#end
#
#
function fill_energy_mesh(band)
    ek = Array{Float32}(undef, 1, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i - 1, j - 1, k - 1, nx, ny, nz)
        ek[1, i, j, k] = band(1, kvec)
    end
    return ek
end
#
#function get_number_of_electrons(G)
#    Ne = 0.0
#    for i in 1:nx, j in 1:ny, k in 1:nz
#        Gsum = sum(real.(G[:, i, j, k]))
#        Ne += 1 + 2 * Gsum / beta
#    end
#    return round(Ne / nk, digits=6)
#end
#
#function calc_electron_density(G, mesh)::Float64
#    """ Calculate electron density from Green function """
#    gio = dropdims(sum(G,dims=(2,3,4)),dims=(2,3,4))/nk
#
#    g_l = fit(mesh.IR_basis_set.smpl_wn_f,gio, dim=1)
#    g_tau0 = dot(mesh.IR_basis_set.basis_f.u(0), g_l)
#
#    n  = 1.0 + real(g_tau0)
#    n  = 2.0 * n #for spin
#    return n
#end
#
#function mu_calc(iw, ek, n, mesh)::Float64
#    """ Find chemical potential for a given filling n0 via brent's root finding algorithm """
#    f = μ -> begin
#        G = 1.0 ./ (iw .- ek .+ μ)  # shape: (fnw, nx, ny, nz)
#        calc_electron_density(G, mesh) - n
#    end
#
#    mu = find_zero(f, (1.1*minimum(ek), 1.1*maximum(ek)), Roots.Brent()) 
#    return mu
#end
#
#
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
#
#
#function get_qvec(ix, iy, iz, nx, ny, nz)
#    kvec = [ix / nx - 0.5, iy / ny - 0.5, iz / nz - 0.5] 
#    if dim < 3
#        kvec[3] = 0.0
#    elseif dim < 2
#        kvec[2] = 0.0
#    end
#    kvec = BZ * kvec
#    return kvec
#end
#
#
#function save(iw, arr, filename)
#    
#    printv("Saving")
#    fnw = length(iw)
#
#    # Determine output header and formatter
#    header, fmt = begin
#        if dim == 3
#            ("# x y z w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f %f %f\n", k[1], k[2], k[3], iw.im, real(val), imag(val)))
#        elseif dim == 2
#            ("# x y w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f %f\n", k[1], k[2], iw.im, real(val), imag(val)))
#        elseif dim == 1
#            ("# x w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f\n", k[1], iw.im, real(val), imag(val)))
#        else
#            error("Invalid dimension: $dim")
#        end
#    end
#
#    # Write data
#    open(filename, "w") do f
#        print(f, header)
#        for ix in 1:nqx, iy in 1:nqy, iz in 1:nqz, l in 1:fnw
#            kvec = get_qvec(ix - 1, iy - 1, iz - 1, nqx, nqy, nqz)
#            w = iw[l]
#            val = arr[l, ix, iy, iz]
#            print(f, fmt(kvec, w, val))
#        end
#    end
#    println("Saved to $filename")
#end



"""
Solver struct to calculate the FLEX loop self-consistently.
After initializing the Solver by `solver = FLEXSolver(mesh, beta, U, n, sigma_init, sfc_tol, maxiter, U_maxiter, mix)'
it can be run by `solve(solver)`.
 """
mutable struct FLEXSolver
    mesh     ::Mesh
    beta     ::Float64
    U        ::Float64
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
end

"""Initiarize function"""
function make_FLEXSolver(
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
        )::FLEXSolver
    
        n::Float64 = 0.0
    
        gkio  = Array{ComplexF32}(undef, mesh.fnw,   nk1, nk2, nk3)
        grit  = Array{ComplexF32}(undef, mesh.fntau, nk1, nk2, nk3)
        ckio  = Array{ComplexF32}(undef, mesh.bnw,   nk1, nk2, nk3)
        V     = Array{ComplexF32}(undef, mesh.bntau, nk1, nk2, nk3)
        sigma = sigma_init
    
        iw, iv = get_iw_iv(mesh)
        iw = reshape(iw, mesh.fnw, 1, 1, 1)
        iv = reshape(iv, mesh.bnw, 1, 1, 1)
        solver = FLEXSolver(mesh, beta, U, n, sfc_tol, maxiter, U_maxiter, mix, verbose, mu, gkio, grit, ckio, V, sigma, iw, iv, ek)
        solver.n = calc_electron_density(solver, mu)
    
        solver.mu = mu_calc(solver)
        gkio_calc(solver,Float32(solver.mu))
        grit_calc(solver)
        ckio_calc(solver)
        return solver
end

#%%%%%%%%%%% Loop solving instance
function solve(solver::FLEXSolver)
    """ FLEXSolver.solve() executes FLEX loop until convergence """
    # check whether U < U_crit! Otherwise, U needs to be renormalized.
    if maximum(abs, solver.ckio) * solver.U >= 1 && scf == true
        U_renormalization(solver)
        println("New U = $(solver.U)")
    end
            
    # perform loop until convergence is reached:
    for it in 1:solver.maxiter
        sigma_old = copy(solver.sigma)
        loop(solver)
        println("max ckio: $(maximum(abs, solver.ckio))")
        println("U: $(solver.U)")
        if maximum(abs, solver.ckio) * solver.U >= 1 
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
    end
end
    
function loop(solver::FLEXSolver)
    """ FLEX loop """
    gkio_old = copy(solver.gkio)
    
    V_calc(solver)
    sigma_calc(solver)
        
    solver.mu = mu_calc(solver)
    gkio_calc(solver,Float32(solver.mu))
    
    solver.gkio .= solver.mix*solver.gkio .+ (1-solver.mix)*gkio_old
        
    grit_calc(solver)
    ckio_calc(solver)
end


#%%%%%%%%%%% U renormalization loop instance
function U_renormalization(solver::FLEXSolver)
    """ Loop for renormalizing U if Stoner enhancement U*max{chi0} >= 1. """
    println("WARNING: U is too large and the spin susceptibility denominator will diverge/turn unphysical!")
    println("Initiate U renormalization loop.")
    
    # save old U for later
    U_old::Float64 = solver.U
    # renormalization loop may run infinitely! Insert break condition after U_it_max steps
    U_it::Int64 = 0
    
    while U_old*maximum(abs, solver.ckio) >= 1.0
        U_it += 1
        # reset U
        solver.U = U_old
        
        # remormalize U such that U*chi0 < 1
        solver.U = solver.U / (maximum(abs, solver.ckio)*solver.U + 0.01)
        println(U_it, '\t', solver.U, '\t', U_old)
        
        # perform one shot FLEX loop
        loop(solver)
        
        
        # break condition for too many steps
        if U_it == solver.U_maxiter
            println("U renormalization reached breaking point")
            break
        end
    end
    println("Leaving U renormalization...")
end

#%%%%%%%%%%% Calculation steps
function gkio_calc(solver::FLEXSolver, mu::Float32)
    solver.gkio .= 1.0 ./ (solver.iw .- (solver.ek .- mu) .- solver.sigma)
    """ calculate Green function G(iw,k) """
end

function grit_calc(solver::FLEXSolver)
    """ Calculate real space Green function G(tau,r) [for calculating chi0 and sigma] """
    solver.grit .= kw_to_rtau(solver.gkio, 'F', solver.mesh)
end

function ckio_calc(solver::FLEXSolver)
    """ Calculate irreducible susciptibility chi0(iv,q) """
    solver.ckio .= rtau_to_kw(solver.grit .* reverse(solver.grit, dims=1), 'B', solver.mesh)
end

function V_FLEX(U, ckio)
    Vkio = (1.5*U^2) .* ckio ./ (1 .- U .* ckio) .+ (0.5*U^2) .* ckio ./ (1 .+ U .* ckio) .- (U^2) .* ckio
    return Vkio
end

function V_calc(solver::FLEXSolver)
    """ Calculate interaction V(tau,r) from RPA-like spin and charge susceptibility for calculating sigma """
    # check whether U is too large and give warning
    maxval = maximum(abs.(solver.ckio))*solver.U
    if maxval >= 1
        error("U*max(chi0) = $(maxval) >= 1! Paramagnetic phase is left and calculations will turn unstable!")
    end

    solver.ckio .= V_FLEX(solver.U, solver.ckio)
    # Constant Hartree Term V ~ U needs to be treated extra, since they cannot be modeled by the IR basis.
    # In the single-band case, the Hartree term can be absorbed into the chemical potential.

    solver.V .= kw_to_rtau(solver.ckio, 'B', solver.mesh)
end

function sigma_calc(solver::FLEXSolver)
    """ Calculate self-energy Sigma(iw,k) """
    solver.grit = solver.V .* solver.grit

    solver.sigma .= rtau_to_kw(solver.grit, 'F', solver.mesh)
end


#%%%%%%%%%%% Setting chemical potential mu
function calc_electron_density(solver::FLEXSolver,mu::Float64)::Float64
    """ Calculate electron density from Green function """
    gkio_calc(solver,Float32(mu))
    gio = dropdims(sum(solver.gkio,dims=(2,3)),dims=(2,3))/solver.mesh.nk

    g_l = fit(solver.mesh.IR_basis_set.smpl_wn_f,gio, dim=1)
    g_tau0 = dot(solver.mesh.IR_basis_set.basis_f.u(0), g_l)

    n  = 1.0 + real(g_tau0)
    n  = 2.0 * n #for spin
end

function mu_calc(solver::FLEXSolver)::Float64
    n = Float64(solver.n)
    """ Find chemical potential for a given filling n0 via brent's root finding algorithm """
    f  = x -> calc_electron_density(solver,Float64(x)) - n

    mu = find_zero(f, (3*Float64(minimum(solver.ek)), 3*Float64(maximum(solver.ek))), Roots.Brent()) 
    return mu
end
# initialize calculation

function main()
    band = Bands()
    ek = fill_energy_mesh(band)
    mesh = IR_Mesh(ek)
    sigma_init = zeros(ComplexF32,(mesh.fnw, nk1, nk2, nk3))
    verbose = cfg.verbosity == "high" 
    solver = make_FLEXSolver(mesh, beta, ek, U, mu, sigma_init, sfc_tol=sfc_tol, maxiter=maxiter, U_maxiter=U_maxiter, mix=mix, verbose=verbose)

    # perform FLEX loop
    solve(solver)
    println("New mu=$(solver.mu)")

    solver.grit = Array{ComplexF32,4}(undef, 0, 0, 0, 0)
    solver.V = Array{ComplexF32,4}(undef, 0, 0, 0, 0)

    V = V_FLEX(solver.U, solver.ckio)

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
end # module
