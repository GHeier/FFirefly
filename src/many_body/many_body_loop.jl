module ManyBodyLoop

using FFTW
using LinearAlgebra
using Roots
using Printf
using Interpolations
using SparseIR
import SparseIR: Statistics, value, valueim
include("../objects/mesh.jl")
using .IRMesh
using Firefly
cfg = Firefly.Config

outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype

scf = cfg.self_consistent 

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


if cfg.interaction != "FLEX"
    println("Code written for FLEX interactions")
    exit()
end

function main()
    band = Bands()
    ek = fill_energy_mesh(band)
    mesh = IR_Mesh(ek)
    iw, iv = get_iw_iv(mesh)
    fnw, bnw, fntau, bntau = mesh.fnw, mesh.bnw, mesh.fntau, mesh.bntau
    iw = reshape(iw, fnw, 1, 1, 1)
    # Preallocate all arrays
    Sigma = zeros(ComplexF32, fnw, nx, ny, nz)
    G     = zeros(ComplexF32, fnw, nx, ny, nz)
    G_rt  = Array{ComplexF32}(undef, fntau, nx, ny, nz)
    X     = Array{ComplexF32}(undef, bnw, nx, ny, nz)
    X_rt  = Array{ComplexF32}(undef, bntau, nx, ny, nz)
    V     = Array{ComplexF32}(undef, bnw, nx, ny, nz)

    loop!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)

    println("Max Self-Energy: $(maximum(abs.(Sigma)))")
    println("Min Self-Energy: $(minimum(abs.(Sigma)))")
    println("Max Vertex: $(maximum(abs.(V)))")
    println("Min Vertex: $(minimum(abs.(V)))")

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
        Sigma[i, :, :, :] .= fftshift(Sigma[i, :, :, :])
    end
    for i in 1:mesh.bnw
        V[i, :, :, :] .= fftshift(V[i, :, :, :])
        X[i, :, :, :] .= fftshift(X[i, :, :, :])
    end

    save_field!(outdir * prefix * "_self_energy." * filetype, Sigma, BZ_in, kmesh, imag.(iw))
    save_field!(outdir * prefix * "_vertex." * filetype, V, BZ_in, kmesh, imag.(iv))
    save_field!(outdir * prefix * "_chi." * filetype, X, BZ_in, kmesh, imag.(iv))
end


function loop!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
    max_iters = 30
    if !scf
        max_iters = 1
    end
    for _ in 1:max_iters
        G .= 1 ./ (iw .- ek .+ mu .- Sigma)

        G_rt .= kw_to_rtau(G, 'F', mesh)

        X_rt .= G_rt .* reverse(G_rt, dims=1)

        X .= rtau_to_kw(X_rt, 'B', mesh)
        max_val = maximum(real.(X)) * U
        if (max_val >= 1)
            println("U*X = $(maximum(abs, X) * U)")
            if scf
                U_renormalization(U, mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
            else
                exit("Diverging Interaction, cannot be handled without scf U renormalization")
            end
        end

        V .= (3/2) .* (U^2 .* X ./ (1 .- U .* X)) .- (1/2) .* (U^2 .* X ./ (1 .+ U .* X))
        #V .= (U^2 .* X ./ (1 .- U .* X)) .+ (U^3 .* X.^2 ./ (1 .+ U^2 .* X.^2))

        X_rt .= kw_to_rtau(V, 'B', mesh)

        G_rt .= X_rt .* G_rt

        Sigma .= rtau_to_kw(G_rt, 'F', mesh)
    end
end

function U_renormalization(U, mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
    """ Loop for renormalizing U if Stoner enhancement U*max{chi0} >= 1. """
    println("WARNING: U is too large and the spin susceptibility denominator will diverge/turn unphysical!")
    println("Initiate U renormalization loop.")
    
    # save old U for later
    U_old::Float64 = U
    # renormalization loop may run infinitely! Insert break condition after U_it_max steps
    U_it::Int64 = 0
    U_maxiter = 50
    
    while U_old*maximum(abs, X) >= 1.0
        U_it += 1
        
        # remormalize U such that U*chi0 < 1
        U = U / (maximum(abs, X)*U + 0.01)
        println(U_it, '\t', U, '\t', U_old)
        
        V .= (3/2) .* (U^2 .* X ./ (1 .- U .* X)) .- (1/2) .* (U^2 .* X ./ (1 .+ U .* X))
        #V .= (U^2 .* X ./ (1 .- U .* X)) .+ (U^3 .* X.^2 ./ (1 .+ U^2 .* X.^2))

        # perform one shot FLEX loop
        loop!(mesh, iw, ek, Sigma, G, G_rt, X, X_rt, V)
        
        # reset U
        U = U_old
        
        # break condition for too many steps
        if U_it == U_maxiter
            println("U renormalization reached breaking point")
            break
        end
    end
    println("Leaving U renormalization...")
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


function get_qvec(ix, iy, iz, nx, ny, nz)
    kvec = [ix / nx - 0.5, iy / ny - 0.5, iz / nz - 0.5] 
    if dim < 3
        kvec[3] = 0.0
    elseif dim < 2
        kvec[2] = 0.0
    end
    kvec = BZ * kvec
    return kvec
end


function fill_energy_mesh(band)
    ek = Array{Float32}(undef, 1, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i - 1, j - 1, k - 1, nx, ny, nz)
        ek[1, i, j, k] = band(1, kvec)
    end
    return ek
end

function calc_electron_density(G, mesh)::Float64
    """ Calculate electron density from Green function """
    gio = dropdims(sum(G,dims=(2,3,4)),dims=(2,3,4))/nk

    g_l = fit(mesh.IR_basis_set.smpl_wn_f,gio, dim=1)
    g_tau0 = dot(mesh.IR_basis_set.basis_f.u(0), g_l)

    n  = 1.0 + real(g_tau0)
    n  = 2.0 * n #for spin
end

function mu_calc(iw, ek, n)::Float64
    """ Find chemical potential for a given filling n0 via brent's root finding algorithm """
    gkio = 1.0 ./ (iw .- ek)
    f  = x -> calc_electron_density(gkio,x) - n

    mu = find_zero(f, (3*minimum(ek), 3*maximum(ek)), Roots.Brent()) 
end

function save(iw, arr, filename)
    
    printv("Saving")
    fnw = length(iw)

    # Determine output header and formatter
    header, fmt = begin
        if dim == 3
            ("# x y z w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f %f %f\n", k[1], k[2], k[3], iw.im, real(val), imag(val)))
        elseif dim == 2
            ("# x y w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f %f\n", k[1], k[2], iw.im, real(val), imag(val)))
        elseif dim == 1
            ("# x w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f\n", k[1], iw.im, real(val), imag(val)))
        else
            error("Invalid dimension: $dim")
        end
    end

    # Write data
    open(filename, "w") do f
        print(f, header)
        for ix in 1:nqx, iy in 1:nqy, iz in 1:nqz, l in 1:fnw
            kvec = get_qvec(ix - 1, iy - 1, iz - 1, nqx, nqy, nqz)
            w = iw[l]
            val = arr[l, ix, iy, iz]
            print(f, fmt(kvec, w, val))
        end
    end
    println("Saved to $filename")
end

end # module
