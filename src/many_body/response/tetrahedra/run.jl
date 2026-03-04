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
    results = Dict{Tuple{Int,Int}, Float64}()

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
    chi_band = zeros(Float64, nw, nqpts)

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
    chi = zeros(Float64, nw, nqpts, nbnd, nbnd, nbnd, nbnd)

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
Only creates non-negative frequencies (w >= 0) to exploit symmetry χ(q, -ω) = χ(q, ω).

# Returns
- Array of non-negative frequency points for calculation
- Total number of frequency points including negative frequencies
"""
function setup_frequency_grid(wpts, w_max)
    if wpts == 1
        return [0.0], 1
    else
        # Calculate number of positive frequencies (including zero if wpts is odd)
        if isodd(wpts)
            # Odd number: include zero and positive frequencies
            n_pos = div(wpts, 2) + 1
            w_pos = collect(range(0.0, w_max, length=n_pos))
        else
            # Even number: only positive frequencies (no zero)
            n_pos = div(wpts, 2)
            dw = 2 * w_max / (wpts - 1)
            w_pos = collect(range(dw, w_max, length=n_pos))
        end
        return w_pos, wpts
    end
end

"""
    mirror_to_negative_frequencies(chi_pos, w_pos, wpts_total)

Mirror response function from positive frequencies to negative frequencies.
Uses symmetry: χ(q, -ω) = χ(q, ω) for real systems.

# Arguments
- `chi_pos`: Response function computed for w >= 0, shape (n_pos, nqpts, nbnd, nbnd, nbnd, nbnd)
- `w_pos`: Positive frequency points
- `wpts_total`: Total number of frequency points (including negative)

# Returns
- Full χ array with negative frequencies filled in, shape (wpts_total, nqpts, ...)
- Full frequency list including negative frequencies
"""
function mirror_to_negative_frequencies(chi_pos, w_pos, wpts_total)
    n_pos = length(w_pos)
    dims = size(chi_pos)
    nqpts = dims[2]
    nbnd_dims = dims[3:end]

    # Create full chi array
    chi_full = zeros(Float64, wpts_total, nqpts, nbnd_dims...)

    # Create full frequency list
    w_max = w_pos[end]
    w_full = collect(range(-w_max, w_max, length=wpts_total))

    if isodd(wpts_total)
        # Odd number of points: w = [..., -dw, 0, dw, ...]
        # Positive frequencies (including zero) go in second half
        i_zero = div(wpts_total, 2) + 1
        chi_full[i_zero:end, :, fill(:, length(nbnd_dims))...] = chi_pos

        # Mirror to negative frequencies (skip zero)
        for i in 1:(i_zero-1)
            i_mirror = wpts_total - i + 1
            chi_full[i, :, fill(:, length(nbnd_dims))...] = chi_pos[i_mirror - i_zero + 1, :, fill(:, length(nbnd_dims))...]
        end
    else
        # Even number of points: w = [..., -dw, dw, ...]
        # Positive frequencies go in second half
        i_mid = div(wpts_total, 2)
        chi_full[(i_mid+1):end, :, fill(:, length(nbnd_dims))...] = chi_pos

        # Mirror to negative frequencies
        for i in 1:i_mid
            i_mirror = wpts_total - i + 1
            chi_full[i, :, fill(:, length(nbnd_dims))...] = chi_pos[i_mirror - i_mid, :, fill(:, length(nbnd_dims))...]
        end
    end

    return chi_full, w_full
end

"""
    find_chi_extrema(chi, w_list, qmesh, BZ, dim; band_indices=(1,1,1,1))

Find the locations of maximum positive and minimum negative values in chi.

# Arguments
- `chi`: Response function array, shape (nw, nqx, nqy, nqz, nbnd, nbnd, nbnd, nbnd)
- `w_list`: Frequency points
- `qmesh`: Q-mesh dimensions (nqx, nqy, nqz)
- `BZ`: Brillouin zone matrix
- `dim`: Spatial dimension (2 or 3)
- `band_indices`: Tuple (i, j, k, l) specifying which band combination to analyze (default: (1,1,1,1))

# Returns
- Dictionary with :max and :min entries, each containing :value, :w, :q, :indices
"""
function find_chi_extrema(chi, w_list, qmesh, BZ, dim; band_indices=(1,1,1,1))
    i, j, k, l = band_indices
    chi_slice = chi[:, :, :, :, i, j, k, l]

    nw = length(w_list)
    nqx, nqy, nqz = qmesh

    # Find max and min values (only consider finite values)
    chi_real = real.(chi_slice)

    # Mask non-finite values
    finite_mask = isfinite.(chi_real)
    if !any(finite_mask)
        @warn "No finite values in chi array"
        return nothing
    end

    # Find maximum positive value
    max_val = -Inf
    max_idx = CartesianIndex(1, 1, 1, 1)
    min_val = Inf
    min_idx = CartesianIndex(1, 1, 1, 1)

    for idx in CartesianIndices(chi_real)
        val = chi_real[idx]
        if isfinite(val)
            if val > max_val
                max_val = val
                max_idx = idx
            end
            if val < min_val
                min_val = val
                min_idx = idx
            end
        end
    end

    # Convert indices to physical coordinates
    function idx_to_coords(idx)
        iw, iqx, iqy, iqz = Tuple(idx)
        w = w_list[iw]

        # Convert q-indices to fractional coordinates (centered grid)
        qfrac = [(iqx - 1) / nqx - 0.5,
                 (iqy - 1) / nqy - 0.5,
                 (iqz - 1) / nqz - 0.5]

        # Convert to Cartesian coordinates using BZ matrix
        q_cart = BZ' * qfrac

        if dim == 2
            q_cart = q_cart[1:2]
            qfrac = qfrac[1:2]
        end

        return (w=w, q_frac=qfrac, q_cart=q_cart, indices=(iw, iqx, iqy, iqz))
    end

    max_coords = idx_to_coords(max_idx)
    min_coords = idx_to_coords(min_idx)

    # Print results
    println("\n" * "="^60)
    println("χ Extrema Analysis (bands: $i,$j,$k,$l)")
    println("="^60)

    println("\n📈 Maximum (most positive):")
    @printf("   Value: %.6f\n", max_val)
    @printf("   Frequency ω: %.4f\n", max_coords.w)
    if dim == 2
        @printf("   q (fractional): (%.4f, %.4f)\n", max_coords.q_frac[1], max_coords.q_frac[2])
        @printf("   q (Cartesian):  (%.4f, %.4f)\n", max_coords.q_cart[1], max_coords.q_cart[2])
    else
        @printf("   q (fractional): (%.4f, %.4f, %.4f)\n", max_coords.q_frac...)
        @printf("   q (Cartesian):  (%.4f, %.4f, %.4f)\n", max_coords.q_cart...)
    end
    @printf("   Indices (iw, iqx, iqy, iqz): %s\n", max_coords.indices)

    println("\n📉 Minimum (most negative):")
    @printf("   Value: %.6f\n", min_val)
    @printf("   Frequency ω: %.4f\n", min_coords.w)
    if dim == 2
        @printf("   q (fractional): (%.4f, %.4f)\n", min_coords.q_frac[1], min_coords.q_frac[2])
        @printf("   q (Cartesian):  (%.4f, %.4f)\n", min_coords.q_cart[1], min_coords.q_cart[2])
    else
        @printf("   q (fractional): (%.4f, %.4f, %.4f)\n", min_coords.q_frac...)
        @printf("   q (Cartesian):  (%.4f, %.4f, %.4f)\n", min_coords.q_cart...)
    end
    @printf("   Indices (iw, iqx, iqy, iqz): %s\n", min_coords.indices)

    println("="^60)

    return Dict(
        :max => (value=max_val, coords=max_coords),
        :min => (value=min_val, coords=min_coords)
    )
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

    # Setup frequency grid (only positive frequencies)
    w_max = cfg.cutoff_energy
    w_pos, wpts_total = setup_frequency_grid(wpts, w_max)
    println("  ω range: [-$(w_max), $(w_max)]")
    println("  Computing for ω ≥ 0 only ($(length(w_pos)) points), will mirror to ω < 0")

    # Load Hamiltonian and compute eigenvalues
    Ek_array, psis, dos, norb = setup_hamiltonian_eigenvalues(kmesh, BZ, dim)

    # Create interpolation grids
    println("\nCreating interpolation grids...")
    Ek_grids, Uk_grids = create_interpolation_grids(Ek_array, psis, kmesh, dim, norb)

    # Calculate response with symmetry reduction (only for w >= 0)
    println("\n" * "="^60)
    println("Starting response grid calculation (ω ≥ 0)")
    println("="^60)
    chi_pos = calculate_response_grid(Ek_grids, Uk_grids, w_pos, kmesh, qmesh, dos, 0)

    # Mirror to negative frequencies using symmetry
    println("\n" * "="^60)
    println("Mirroring to negative frequencies using χ(q,-ω) = χ(q,ω)")
    println("="^60)
    chi, w_full = mirror_to_negative_frequencies(chi_pos, w_pos, wpts_total)
    println("✓ Full frequency grid constructed: $(wpts_total) points")

    # Analyze extrema
    extrema_info = find_chi_extrema(chi, w_full, qmesh, BZ, dim)

    # Save results
    save_response_results(chi, w_full, qmesh, BZ, dim)

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
