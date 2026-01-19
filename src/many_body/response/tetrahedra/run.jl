using Firefly
cfg = Firefly.Config

# Load relevant variables from the configuration

using LinearAlgebra
using Printf
using Interpolations
using Interpolations: extrapolate, Periodic
using Base.Threads

# Import utility functions
include("response_utils.jl")

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
if cfg.mu_from_n
    n_field = Field_R(outdir*prefix*"_E_vs_n."*filetype)
    mu = n_field(cfg.num_electrons)
    println("Shifting mu to $(mu) based on electron number $(cfg.num_electrons)")
end
U = cfg.U0
BZ = cfg.brillouin_zone
celltype = cfg.celltype

"""
    setup_symmetry_reduction(qmesh, celltype, dim)

Set up symmetry reduction for q-point grid.

# Arguments
- `qmesh`: Q-mesh dimensions (nqx, nqy, nqz)
- `celltype`: Cell/lattice type (e.g., "SC", "BCC", "FCC")
- `dim`: Spatial dimension (2 or 3)

# Returns
- `unique_to_qindices`: Vector mapping group index → list of equivalent q-point indices
- `nqpts_unique`: Number of unique q-point groups
- `nqpts`: Total number of q-points
"""
function setup_symmetry_reduction(qmesh, celltype, dim)
    nqpts = prod(qmesh)

    println("\nApplying symmetry reduction to q-grid...")
    reduced_qgrid = get_reduced_grid(collect(qmesh), celltype)
    nqpts_unique = length(reduced_qgrid)
    reduction_pct = round(100 * (1 - nqpts_unique / nqpts), digits=1)

    println("  Cell type: $(celltype)")
    println("  Unique q-points: $nqpts_unique / $nqpts ($reduction_pct% reduction)")

    # Create mapping: q-index → unique group index
    q_to_group = zeros(Int, nqpts)
    unique_to_qindices = Vector{Vector{Int}}(undef, nqpts_unique)

    for (igroup, group) in enumerate(reduced_qgrid)
        unique_to_qindices[igroup] = []
        for point in group
            # Convert [ix, iy, iz] indices to linear index (0-based to 1-based)
            iq = 1 + point[1] + point[2] * qmesh[1]
            if dim == 3
                iq += point[3] * qmesh[1] * qmesh[2]
            end
            q_to_group[iq] = igroup
            push!(unique_to_qindices[igroup], iq)
        end
    end

    # Verify mapping
    if length(unique(q_to_group)) != nqpts_unique
        error("Symmetry mapping error: not all unique groups represented")
    end
    println("  ✓ Symmetry mapping verified")

    return unique_to_qindices, nqpts_unique, nqpts
end

"""
    evaluate_grid_for_dimension(ipt, kmesh, dim)

Evaluate interpolation object on grid for given dimension.

# Arguments
- `ipt`: Interpolation object
- `kmesh`: K-mesh dimensions (nkx, nky, nkz)
- `dim`: Spatial dimension (2 or 3)

# Returns
- Grid evaluation result
"""
function evaluate_grid_for_dimension(ipt, kmesh, dim)
    nkx, nky, nkz = kmesh
    if dim == 2
        return ipt(1:nkx, 1:nky)
    else
        return ipt(1:nkx, 1:nky, 1:nkz)
    end
end

"""
    calculate_qw_response(igroup, unique_to_qindices, qpts, w_pts, Ek_grid, Ek_grids,
                         k, l, dos, mu, kmesh, dim, iter)

Calculate response for all ω at a single q-group.

# Returns
- Dictionary mapping (iw, iq) → response value
"""
function calculate_qw_response(igroup, unique_to_qindices, qpts, w_pts,
                               Ek_grid, Ek_grids, k, l, dos, mu, kmesh, dim, iter)
    results = Dict{Tuple{Int,Int}, ComplexF64}()

    # Get representative q-point
    iq_rep = unique_to_qindices[igroup][1]
    q = qpts[iq_rep, :]

    # Evaluate E(k+q) once for this q-group
    Ekq_grid = evaluate_Ekq_on_grid(Ek_grids[k, l], kmesh, q, dim)

    for (iw, w) in enumerate(w_pts)
        # Special case: q=0, w=0 → static susceptibility = DOS
        if w == 0.0
            q_norm = q[1]^2 + q[2]^2 + (dim == 3 ? q[3]^2 : 0)
            if q_norm < 1e-6
                for iq in unique_to_qindices[igroup]
                    results[(iw, iq)] = dos
                end
                continue
            end
        end

        # Calculate response for representative q-point
        result = calculate_response_bzintegral_2(w, Ek_grid, Ekq_grid, mu, iter)

        # Store for all symmetry-equivalent q-points
        for iq in unique_to_qindices[igroup]
            results[(iw, iq)] = result
        end
    end

    return results
end

"""
    calculate_band_pair_response(i, j, k, l, Ek_grids, w_pts, unique_to_qindices,
                                 qpts, nqpts_unique, dos, mu, kmesh, dim, iter,
                                 progress_counter, progress_lock, total_iterations, start_time)

Calculate response for one band pair combination.

# Returns
- 4D array chi[iw, iq] for the given band indices
"""
function calculate_band_pair_response(i, j, k, l, Ek_grids, w_pts, unique_to_qindices,
                                     qpts, nqpts_unique, nw, nqpts, dos, mu, kmesh, dim, iter,
                                     progress_counter, progress_lock, total_iterations, start_time)
    # Evaluate E_ij(k) once
    Ek_grid = evaluate_grid_for_dimension(Ek_grids[i, j], kmesh, dim)

    # Allocate result array for this band pair
    chi_band = zeros(ComplexF64, nw, nqpts)

    # Parallelize over unique q-groups
    Threads.@threads for igroup in 1:nqpts_unique
        results = calculate_qw_response(igroup, unique_to_qindices, qpts, w_pts,
                                       Ek_grid, Ek_grids, k, l, dos, mu, kmesh, dim, iter)

        # Store results
        for ((iw, iq), value) in results
            chi_band[iw, iq] = value
        end

        # Thread-safe progress update
        current = atomic_add!(progress_counter, length(w_pts))
        if current % 10 == 0
            lock(progress_lock) do
                print_progress(current, total_iterations, start_time)
            end
        end
    end

    return chi_band
end

"""
    calculate_response_grid(Ek_grids, Uk_grids, w_pts, kmesh, qmesh, dos, iter=2)

Calculate response function on q-ω grid with symmetry reduction.

# Returns
- 8D array chi[iw, iqx, iqy, iqz, i, j, k, l] with response function
"""
function calculate_response_grid(Ek_grids, Uk_grids, w_pts, kmesh, qmesh, dos, iter=2)
    nw = length(w_pts)
    nkx, nky, nkz = kmesh
    nqx, nqy, nqz = qmesh
    qpts = get_fractional_mesh(qmesh; centered=true)

    nbnd = Int(sqrt(length(Ek_grids)))
    println("nbnd: ", nbnd)

    # Setup symmetry reduction
    unique_to_qindices, nqpts_unique, nqpts = setup_symmetry_reduction(qmesh, celltype, dim)

    # Initialize susceptibility
    chi = zeros(ComplexF64, nw, nqpts, nbnd, nbnd, nbnd, nbnd)

    # Setup progress tracking
    total_iterations = nbnd^4 * nqpts_unique * nw
    progress_counter = Atomic{Int}(0)
    progress_lock = ReentrantLock()

    println("\nStarting band and q-ω loops...")
    println("Using ", nthreads(), " threads for parallel computation")
    println("Total iterations: $total_iterations (reduced from $(nbnd^4 * nqpts * nw))")
    start_time = time()

    # Loop over all band combinations
    for i in 1:nbnd, j in 1:nbnd, k in 1:nbnd, l in 1:nbnd
        chi_band = calculate_band_pair_response(i, j, k, l, Ek_grids, w_pts,
                                                unique_to_qindices, qpts, nqpts_unique,
                                                nw, nqpts, dos, mu, kmesh, dim, iter,
                                                progress_counter, progress_lock,
                                                total_iterations, start_time)

        # Store in full chi array
        chi[:, :, i, j, k, l] = chi_band

        println("\nCompleted bands (", i, ", ", j, ", ", k, ", ", l, ")")
    end

    # Print summary
    println("\n✓ Grid calculation complete!")
    elapsed_total = time() - start_time
    speedup = (nbnd^4 * nqpts * nw) / total_iterations
    println("Total computation time: $(round(elapsed_total, digits=3)) s")
    println("Effective speedup from symmetry: $(round(speedup, digits=2))x")
    println()

    # Reshape to 8D
    chi = reshape(chi, (nw, nqx, nqy, nqz, nbnd, nbnd, nbnd, nbnd))
    return chi
end

"""
    setup_hamiltonian_eigenvalues(kmesh, BZ, dim)

Load Hamiltonian and compute eigenvalues/wavefunctions on k-mesh.

# Returns
- `Ek_array`: Eigenvalues array [nkpts, norb, norb]
- `psis`: Wavefunctions array [nkpts, norb, norb]
- `dos`: Density of states at Fermi level
- `norb`: Number of orbitals/bands
"""
function setup_hamiltonian_eigenvalues(kmesh, BZ, dim)
    println("Loading Hamiltonian...")
    kpts = get_kmesh(BZ, kmesh; centered=true)
    H = Firefly.Hamiltonian()

    println("Computing eigenvalues and wavefunctions on k-mesh...")
    Ek_array, psis = get_wavefunctions(H, kpts)

    # Calculate DOS at Fermi level
    Ek_mesh = reshape(real(Ek_array[:, 1, 1]), kmesh...)
    if dim == 2
        Ek_mesh = dropdims(Ek_mesh, dims=3)
    end
    dos = calculate_dos(mu, Ek_mesh, 0)

    norb = size(Ek_array, 2)
    nkpts = size(Ek_array, 1)

    println("  DOS at μ = $(mu): $(dos) states/unit cell/eV")
    println("  nbnd: $(norb), nkpts: $(nkpts)")

    return Ek_array, psis, dos, norb
end

"""
    create_interpolation_grids(Ek_array, psis, kmesh, dim, norb)

Create periodic interpolation grids for eigenvalues and wavefunctions.

# Returns
- `Ek_grids`: Matrix of interpolation objects for E(k)
- `Uk_grids`: Matrix of interpolation objects for wavefunctions
"""
function create_interpolation_grids(Ek_array, psis, kmesh, dim, norb)
    nkx, nky, nkz = kmesh

    Ek_grids = Matrix{Interpolations.AbstractInterpolation}(undef, norb, norb)
    Uk_grids = Matrix{Interpolations.AbstractInterpolation}(undef, norb, norb)

    for i in 1:norb, j in 1:norb
        Ek_band = reshape(real(Ek_array[:, i, j]), nkx, nky, nkz)
        Uk = reshape(psis[:, i, j], nkx, nky, nkz)

        if dim == 2
            Ek_band = dropdims(Ek_band, dims=3)
            Uk = dropdims(Uk, dims=3)
        end

        # Create interpolation with periodic boundary conditions
        itp_Ek = interpolate(Ek_band, BSpline(Linear()))
        itp_Uk = interpolate(Uk, BSpline(Linear()))
        itp_Ek_periodic = extrapolate(itp_Ek, Periodic())
        itp_Uk_periodic = extrapolate(itp_Uk, Periodic())

        println("Created interpolation for bands (", i, ", ", j, ")")

        Ek_grids[i, j] = itp_Ek_periodic
        Uk_grids[i, j] = itp_Uk_periodic
    end

    return Ek_grids, Uk_grids
end

"""
    setup_frequency_grid(wpts, w_max)

Create frequency grid for response calculation.

# Returns
- Array of frequency points
"""
function setup_frequency_grid(wpts, w_max)
    if wpts == 1
        return [0.0]
    else
        return collect(range(-w_max, w_max, length=wpts))
    end
end

"""
    save_response_results(chi, w_list, qmesh, BZ, dim)

Save response function results to file.
"""
function save_response_results(chi, w_list, qmesh, BZ, dim)
    println("\nResponse function statistics:")
    println("  Max real part: ", maximum(x -> isfinite(x) ? x : -Inf, real(chi)))
    println("  Min real part: ", minimum(x -> isfinite(x) ? x : Inf, real(chi)))
    println("  Max chi(w_max): ", maximum(real(chi[end, :, :, :, 1, 1, 1, 1])))
    println("  Min chi(w_max): ", minimum(real(chi[end, :, :, :, 1, 1, 1, 1])))

    BZ_save = BZ[1:dim, 1:dim]
    qmesh_save = qmesh[1:dim]

    filename = outdir*prefix*"_chi."*filetype
    save_data!(filename, chi[:,:,:,:,1,1,1,1], qmesh_save, BZ_save, w_points=w_list)
    println("Saved response function to ", filename)
end

"""
    response_bz_integral()

Main driver function for response function calculation.

Orchestrates the entire calculation:
1. Sets up frequency grid
2. Loads Hamiltonian and computes eigenvalues
3. Creates interpolation grids
4. Calculates response with symmetry reduction
5. Saves results

# Returns
- Maximum finite real part of susceptibility
"""
function response_bz_integral()
    println("="^60)
    println("Calculating response function using BZIntegral")
    println("="^60)
    println("  nbnd: $(nbnd)")
    println("  μ: $(mu)")
    println("  U: $(U)")
    println("  Dimension: $(dim)D")
    println("  w-pts: $(wpts)")
    println("  k-mesh: $(kmesh)")
    println("  q-mesh: $(qmesh)")

    # Setup frequency grid
    w_max = cfg.cutoff_energy
    w_list = setup_frequency_grid(wpts, w_max)
    println("  ω range: [-$(w_max), $(w_max)]")

    # Load Hamiltonian and compute eigenvalues
    Ek_array, psis, dos, norb = setup_hamiltonian_eigenvalues(kmesh, BZ, dim)

    # Create interpolation grids
    println("\nCreating interpolation grids...")
    Ek_grids, Uk_grids = create_interpolation_grids(Ek_array, psis, kmesh, dim, norb)

    # Calculate response with symmetry reduction
    println("\n" * "="^60)
    println("Starting response grid calculation")
    println("="^60)
    chi = calculate_response_grid(Ek_grids, Uk_grids, w_list, kmesh, qmesh, dos, 0)

    # Save results
    save_response_results(chi, w_list, qmesh, BZ, dim)

    return maximum(x -> isfinite(x) ? x : -Inf, real(chi))
end

"""
    run()

Entry point for response function calculation.
"""
function run()
    chi_max = response_bz_integral()
    return chi_max
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
