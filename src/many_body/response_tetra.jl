#!/usr/bin/env julia
#
# Using BZIntegral.jl (Tetrahedron Method)
#
module Tetrahedron

using LinearAlgebra
using Printf
using BZIntegral
using BZIntegral.BZInt2D
using Interpolations
using Interpolations: extrapolate, Periodic
using Base.Threads

using Firefly
cfg = Firefly.Config

outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype

kmesh = cfg.k_mesh
# Ensure odd mesh size for BZIntegral
for i in 1:length(kmesh)
    if iseven(kmesh[i])
        kmesh[i] += 1
    end
end

qmesh = cfg.q_mesh
dim = cfg.dimension
if dim == 2
    kmesh[3] = 1
    qmesh[3] = 1
end

wpts = cfg.w_pts
nbnd = cfg.nbnd
mu = cfg.fermi_energy
U = cfg.onsite_U
BZ = cfg.brillouin_zone


function get_fractional_mesh(n; centered=true)
    nx, ny, nz = n

    if centered
        xs = (0:nx-1) ./ nx .- 0.5
        ys = (0:ny-1) ./ ny .- 0.5
        zs = (0:nz-1) ./ nz .- 0.5
    else
        xs = (0:nx-1) ./ nx
        ys = (0:ny-1) ./ ny
        zs = (0:nz-1) ./ nz
    end

    # fractional k-points: 3 × Nk
    frac = hcat(vec(repeat(xs, inner=(ny*nz))),
                vec(repeat(ys, inner=(nz,), outer=(nx))),
                vec(repeat(zs, outer=(nx*ny))))'
    # Return as Nk × 3 matrix (transpose from 3 × Nk)
    return Matrix{Float64}(frac')
end

function get_kmesh(BZ::AbstractMatrix{<:Real}, n; centered=true)
    nx, ny, nz = n

    if centered
        xs = (0:nx-1) ./ nx .- 0.5
        ys = (0:ny-1) ./ ny .- 0.5
        zs = (0:nz-1) ./ nz .- 0.5
    else
        xs = (0:nx-1) ./ nx
        ys = (0:ny-1) ./ ny
        zs = (0:nz-1) ./ nz
    end

    # fractional k-points: 3 × Nk
    frac = hcat(vec(repeat(xs, inner=(ny*nz))),
                vec(repeat(ys, inner=(nz,), outer=(nx))),
                vec(repeat(zs, outer=(nx*ny))))'
    # Cartesian k-points: 3 × Nk
    kpts = BZ * frac
    # Return as Nk × 3 matrix (transpose from 3 × Nk)
    return Matrix{Float64}(kpts')
end

function evaluate_Ekq_on_grid(itp_Ek, kmesh, q_frac, dim=3)
    """
    Evaluate E(k+q) for all k on the grid with periodic boundary conditions.

    Args:
        itp_Ek: Interpolation object from Interpolations.jl
        kmesh: (nkx, nky, nkz) grid dimensions
        q_frac: [qx, qy, qz] momentum transfer in fractional coordinates
        dim: spatial dimension (2 or 3)

    Returns:
        Array of E(k+q) values on the k-grid with shape (nkx, nky, nkz) or (nkx, nky)
    """
    nkx, nky, nkz = kmesh

    # Original k-grid in fractional coordinates (centered: -0.5 to 0.5)
    kx_frac = (0:nkx-1) ./ nkx .- 0.5
    ky_frac = (0:nky-1) ./ nky .- 0.5
    kz_frac = (0:nkz-1) ./ nkz .- 0.5

    # Add q shift and apply periodic boundary conditions: wrap to [-0.5, 0.5)
    kxq_frac = mod.(kx_frac .+ q_frac[1] .+ 0.5, 1.0) .- 0.5
    kyq_frac = mod.(ky_frac .+ q_frac[2] .+ 0.5, 1.0) .- 0.5
    kzq_frac = mod.(kz_frac .+ q_frac[3] .+ 0.5, 1.0) .- 0.5

    # Convert to interpolation indices (1-based)
    Ikxq = (kxq_frac .+ 0.5) .* nkx .+ 1
    Ikyq = (kyq_frac .+ 0.5) .* nky .+ 1
    Ikzq = (kzq_frac .+ 0.5) .* nkz .+ 1

    # Evaluate interpolation with broadcasting
    if dim == 2
        Ix = repeat(reshape(Ikxq, :, 1), 1, nky)
        Iy = repeat(reshape(Ikyq, 1, :), nkx, 1)
        return itp_Ek.(Ix, Iy)
    else
        Ix = repeat(reshape(Ikxq, :, 1, 1), 1, nky, nkz)
        Iy = repeat(reshape(Ikyq, 1, :, 1), nkx, 1, nkz)
        Iz = repeat(reshape(Ikzq, 1, 1, :), nkx, nky, 1)
        return itp_Ek.(Ix, Iy, Iz)
    end
end

"""
    calculate_response_bzintegral_2(qx, qy, w, Ek_mesh, Ekq_mesh, kmesh, beta, eta, vol)

Calculate response function using BZIntegral tetrahedron method.

χ₀(q,ω) = -2 ∫dk [f(ε(k+q)) - f(ε(k))] / [ω + ε(k) - ε(k+q) + iη]

For the integrand F(k) = 1/[ω + ε(k) - ε(k+q) + iη], we compute:
  ∫ Θ(-ε(k+q)) * F(k) dk  -  ∫ Θ(-ε(k)) * F(k) dk
"""
function calculate_response_bzintegral_2(w, Ek_mesh, Ekq_mesh, mu, iter=2)
    # Denominator: ω + ε(k) - ε(k+q) + small eta to avoid division by zero
    # Using real eta since Quad2DRuleΘ𝔇 doesn't support complex denominators
    eta = 1e-4  # Small broadening parameter
    denom = w .+ Ek_mesh .- Ekq_mesh .+ eta

    Wmesh = Quad2DRuleΘ𝔇(Ek_mesh, mu , denom, iter)-Quad2DRuleΘ𝔇(Ekq_mesh, mu, denom, iter)
    out = -2.0 * sum(Wmesh)  # Include -2 spin factor
    return out
end

function print_progress(current, total, start_time)
    elapsed = time() - start_time
    rate = current / elapsed
    remaining = (total - current) / rate
    @printf("  Progress: %d/%d | Rate: %.1f/s | ETA: %.1f s\r", current, total, rate, remaining)
    flush(stdout)
end

"""
    calculate_response_grid(kmesh, qmesh, w_pts, iter=2)

Efficiently calculate response function on a q-ω grid.
Pre-computes Ek_mesh, then loops over q-points and ω-points.

Returns: 3D array chi[iw, iqx, iqy]
"""
function calculate_response_grid(Ek_grids, Uk_grids, w_pts, kmesh, qmesh, BZ, iter=2)
    nw = length(w_pts)
    nkx, nky, nkz = kmesh
    nqx, nqy, nqz = qmesh
    nqpts = prod(qmesh)
    qpts = get_fractional_mesh(qmesh; centered=true)

    # Get nbnd from Ek_grids
    nbnd = Int(sqrt(length(Ek_grids)))
    println("nbnd: ", nbnd)

    # Initialize susceptibility
    chi = zeros(ComplexF64, nw, nqpts, nbnd, nbnd, nbnd, nbnd)
    println("Starting band and q-ω loops...")
    start_time = time()

    function eval(ipt)
        if dim == 2
            return ipt(1:nkx, 1:nky)
        else
            return ipt(1:nkx, 1:nky, 1:nkz)
        end
    end

    # Total iterations for progress tracking
    total_iterations = nbnd^4 * nqpts * nw
    progress_counter = Atomic{Int}(0)
    progress_lock = ReentrantLock()

    println("Using ", nthreads(), " threads for parallel computation")

    for i in 1:nbnd, j in 1:nbnd
        Ek_grid = eval(Ek_grids[i, j])
        for k in 1:nbnd, l in 1:nbnd
            # Parallelize over q-points
            Threads.@threads for iq in 1:nqpts
            #for iq in 1:nqpts
                q = qpts[iq, :]
                #println("q = ", q)
                Ekq_grid = evaluate_Ekq_on_grid(Ek_grids[k, l], (nkx, nky, nkz), q, dim)
                for (iw, w) in enumerate(w_pts)
                    result = calculate_response_bzintegral_2(w, Ek_grid, Ekq_grid, mu, iter)
                    chi[iw, iq, i, j, k, l] = result

                    # Thread-safe progress update
                    current = atomic_add!(progress_counter, 1)

                    # Only print progress from one thread occasionally
                    if current % 10 == 0
                        lock(progress_lock) do
                            print_progress(current, total_iterations, start_time)
                        end
                    end
                end
            end
            println("\nCompleted bands (", i, ", ", j, ", ", k, ", ", l, ")")
        end
    end

    println("\n✓ Grid calculation complete!")
    elapsed_total = time() - start_time
    println("Total computation time: $(round(elapsed_total, digits=3)) s")
    println()
    chi = reshape(chi, (nw, nqx, nqy, nqz, nbnd, nbnd, nbnd, nbnd))

    return chi
end

function response_bz_integral()
    print("Calculating response function using BZIntegral (Tetrahedron Method)...\n")
    print("  nbnd: $(nbnd)\n")
    print("  μ: $(mu)\n")
    print("  U: $(U)\n")
    print("  Dimension: $(dim)D\n")
    print("  w-pts: $(wpts)\n")
    w_max = cfg.cutoff_energy
    if wpts == 1
        w_list = [0.0]
    else
        w_list = collect(range(-w_max, w_max, length=wpts))
    end
    #w_list = range(-w_max, w_max, length=wpts)
    print("  Taking ω ∈ [0, $(w_max)] and projecting for -ω\n")

    println("Setting up grid calculation:")
    print("  k-mesh: $(kmesh)\n")
    print("  q-mesh: $(qmesh)\n")
    println()

    # Create k-mesh and evaluate Hamiltonian once
    println("Loading Hamiltonian...")
    kpts = get_kmesh(BZ, kmesh; centered=true)
    H = Firefly.Hamiltonian()
    println("Computing eigenvalues and wavefunctions on k-mesh...")
    Ek_array, psis = get_wavefunctions(H, kpts)
    #kpts = get_fractional_mesh(kmesh; centered=true)
    #qpts = get_fractional_mesh(qmesh; centered=true)
    qpts = get_kmesh(BZ, qmesh; centered=true)

    println("Eig shape = ", size(Ek_array))
    println("Psi shape = ", size(psis))
    nkpts = size(Ek_array, 1)
    norb = size(Ek_array, 2)
    if nbnd != norb
        println("  Warning: nbnd from config ($(nbnd)) does not match Hamiltonian output ($(norb)). Using $(norb).")
    end
    nkx, nky, nkz = kmesh
    println("  nbnd: $(norb), nkpts: $(nkpts)")

    Ek_grids = Matrix{Interpolations.AbstractInterpolation}(undef, norb, norb)
    Uk_grids = Matrix{Interpolations.AbstractInterpolation}(undef, norb, norb)

    for i in 1:norb, j in 1:norb
        Ek_band = reshape(real(Ek_array[:, i, j]), nkx, nky, nkz)
        Uk = reshape(psis[:, i, j], nkx, nky, nkz)
        if dim == 2
            Ek_band = dropdims(Ek_band, dims=3)  # Remove z-dimension for 2D
            Uk = dropdims(Uk, dims=3)
        end

        itp_Ek = interpolate(Ek_band, BSpline(Linear()))
        itp_Uk = interpolate(Uk, BSpline(Linear()))

        # Wrap with periodic boundary conditions
        itp_Ek_periodic = extrapolate(itp_Ek, Periodic())
        itp_Uk_periodic = extrapolate(itp_Uk, Periodic())

        println("Created interpolation for bands (", i, ", ", j, ")")

        Ek_grids[i, j] = itp_Ek_periodic
        Uk_grids[i, j] = itp_Uk_periodic
    end

    println("\nStarting response grid calculation...")
    # Call calculate_response_grid with energy grids
    chi = calculate_response_grid(Ek_grids, Uk_grids, w_list, kmesh, qmesh, BZ, 2)
    #print(chi)
    println("Max chi real part: ", maximum(x -> isfinite(x) ? x : -Inf, real(chi)))
    println("Min chi real part: ", minimum(x -> isfinite(x) ? x : -Inf, real(chi)))
    println("Max chi(w_max) = ", maximum(real(chi[end, :, :, :, 1, 1, 1, 1])))
    println("Min chi(w_max) = ", minimum(real(chi[end, :, :, :, 1, 1, 1, 1])))
    chi_neg = conj.(chi)
    #chi = vcat(reverse(chi_neg[2:end, :], dims=1), chi)  # Combine -ω and +ω
    println("w_points = ", w_list)
    BZ_save = BZ[1:dim, 1:dim]  # Adjust BZ for dimension
    qmesh_save = qmesh[1:dim]  # Adjust qmesh for dimension
    save_data!(outdir*prefix*"_chi."*filetype, chi[:,:,:,:,1,1,1,1], qmesh_save, BZ_save, w_points=w_list)
    println("Saved response function to "*outdir*prefix*"_chi."*filetype)
    
end


end # module
