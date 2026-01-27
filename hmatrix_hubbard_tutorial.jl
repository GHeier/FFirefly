#!/usr/bin/env julia
# =============================================================================
# Simple HMatrix Tutorial: Hubbard Model Hamiltonian Diagonalization
# =============================================================================
#
# This demonstrates:
# 1. Building the Hubbard Hamiltonian matrix in Fock space (many-body basis)
# 2. Using sparse matrices for efficient storage
# 3. Finding ground state energy with Arnoldi iteration
#
# Physics: Solve H|ψ⟩ = E|ψ⟩ for the Hubbard model
# H = -t Σ_⟨i,j⟩,σ (c†_{i,σ} c_{j,σ} + h.c.) + U Σ_i n_{i,↑} n_{i,↓}
#
# where:
#   - First term: kinetic energy (hopping between nearest neighbors)
#   - Second term: on-site Coulomb repulsion (Hubbard U)
#   - σ = ↑,↓ (spin)
#
# The Hilbert space dimension grows as C(N_sites, N_up) × C(N_sites, N_down).
# For half-filling (N_up = N_down = N_sites/2), this grows very quickly!
#
# Examples:
#   4×2 = 8 sites:  C(8,4)² = 70² = 4,900
#   4×3 = 12 sites: C(12,6)² = 924² = 853,776
#   4×4 = 16 sites: C(16,8)² = 12,870² = 165,635,700
#   6×6 = 36 sites: C(36,18)² ≈ 8.2 × 10¹⁹ (too large!)
#
# For systems beyond ~16 sites, consider using symmetry reduction (momentum,
# spin, particle-hole) or other methods like DMRG or quantum Monte Carlo.

using HMatrices
using KrylovKit
using LinearAlgebra
using SparseArrays
using StaticArrays
using Printf

# =============================================================================
# PARAMETERS
# =============================================================================

# Physical parameters
const t = 1.0           # Hopping parameter (eV)
const U = 4.0           # Hubbard U (eV)

# System size (start small for tutorial)
const Lx = 4            # Lattice size in x
const Ly = 2            # Lattice size in y
const N_sites = Lx * Ly # Total number of sites (8 sites)
                        # Increase Ly to 3 for 12 sites (dim = 853,776)
                        # or Ly = 4 for 16 sites (dim = 165,635,700)

# Number of particles
const N_up = N_sites ÷ 2      # Spin-up electrons (half-filling)
const N_down = N_sites ÷ 2    # Spin-down electrons

# Numerical parameters
const n_eigenvalues = 3       # Number of lowest eigenvalues to compute

println("="^70)
println("HMatrix Tutorial: Hubbard Model Hamiltonian")
println("="^70)
println("\nSystem parameters:")
@printf("  Lattice: %d × %d (N = %d sites)\n", Lx, Ly, N_sites)
@printf("  Hopping: t = %.2f eV\n", t)
@printf("  Hubbard U: U = %.2f eV\n", U)
@printf("  N_up = %d, N_down = %d (half-filling)\n", N_up, N_down)

# =============================================================================
# 1. BASIS CONSTRUCTION
# =============================================================================

"""
Represent a many-body state as a pair of integers (basis_up, basis_down)
where each bit represents occupation of a site.
For example, if bit i is set, site i is occupied.
"""
struct ManyBodyBasis
    states_up::Vector{Int}      # All possible spin-up configurations
    states_down::Vector{Int}    # All possible spin-down configurations
    n_up::Int
    n_down::Int
    n_sites::Int
end

"""
Generate all N-particle states on n_sites (using bit representation).
This is equivalent to choosing N sites out of n_sites: C(n_sites, N).
"""
function generate_fock_states(n_sites::Int, n_particles::Int)
    states = Int[]

    function generate_recursive(current_state, site, particles_left)
        if particles_left == 0
            push!(states, current_state)
            return
        end
        if site > n_sites || (n_sites - site + 1) < particles_left
            return
        end

        # Don't put particle at this site
        generate_recursive(current_state, site + 1, particles_left)

        # Put particle at this site
        new_state = current_state | (1 << (site - 1))
        generate_recursive(new_state, site + 1, particles_left - 1)
    end

    generate_recursive(0, 1, n_particles)
    return states
end

"""
Create the many-body basis for the Hubbard model.
"""
function create_basis(n_sites, n_up, n_down)
    println("\n" * "-"^70)
    println("Building Fock space basis...")
    println("-"^70)

    states_up = generate_fock_states(n_sites, n_up)
    states_down = generate_fock_states(n_sites, n_down)

    dim_up = length(states_up)
    dim_down = length(states_down)
    dim_total = dim_up * dim_down

    @printf("  Spin-up Hilbert space dimension: %d\n", dim_up)
    @printf("  Spin-down Hilbert space dimension: %d\n", dim_down)
    @printf("  Total Hilbert space dimension: %d\n", dim_total)

    # Calculate expected dimension: C(N_sites, N_up) × C(N_sites, N_down)
    expected = binomial(n_sites, n_up) * binomial(n_sites, n_down)
    @printf("  Expected: %d (match: %s)\n", expected, dim_total == expected)

    return ManyBodyBasis(states_up, states_down, n_up, n_down, n_sites)
end

# =============================================================================
# 2. HAMILTONIAN MATRIX ELEMENTS
# =============================================================================

"""
Check if site i is occupied in state (returns 0 or 1).
"""
@inline function is_occupied(state::Int, site::Int)
    return (state >> (site - 1)) & 1
end

"""
Create electron at site i (returns new state and sign, or nothing if already occupied).
"""
@inline function create_electron(state::Int, site::Int)
    if is_occupied(state, site) == 1
        return nothing  # Site already occupied
    end

    # Count number of particles to the left (for fermionic sign)
    mask = (1 << (site - 1)) - 1
    n_left = count_ones(state & mask)
    sign = (-1)^n_left

    new_state = state | (1 << (site - 1))
    return (new_state, sign)
end

"""
Annihilate electron at site i.
"""
@inline function annihilate_electron(state::Int, site::Int)
    if is_occupied(state, site) == 0
        return nothing  # Site not occupied
    end

    mask = (1 << (site - 1)) - 1
    n_left = count_ones(state & mask)
    sign = (-1)^n_left

    new_state = state & ~(1 << (site - 1))
    return (new_state, sign)
end

"""
Apply hopping operator c†_j c_i (hop from site i to site j) with spin σ.
Returns (new_state, matrix_element) or nothing.
"""
function apply_hopping(state::Int, site_i::Int, site_j::Int)
    # Annihilate at i
    result_i = annihilate_electron(state, site_i)
    if result_i === nothing
        return nothing
    end
    intermediate_state, sign_i = result_i

    # Create at j
    result_j = create_electron(intermediate_state, site_j)
    if result_j === nothing
        return nothing
    end
    final_state, sign_j = result_j

    return (final_state, sign_i * sign_j)
end

"""
Get nearest neighbors on 2D lattice with periodic boundary conditions.
Site indexing: site = x + y * Lx, where x ∈ [0, Lx-1], y ∈ [0, Ly-1].
"""
function get_neighbors(site::Int, Lx::Int, Ly::Int)
    # Convert to (x, y) coordinates (0-indexed)
    x = (site - 1) % Lx
    y = (site - 1) ÷ Lx

    # Neighbors with periodic boundary conditions
    right = ((x + 1) % Lx) + y * Lx + 1
    left = ((x - 1 + Lx) % Lx) + y * Lx + 1
    up = x + ((y + 1) % Ly) * Lx + 1
    down = x + ((y - 1 + Ly) % Ly) * Lx + 1

    return [right, left, up, down]
end

"""
Compute Hamiltonian matrix element ⟨ψ_i|H|ψ_j⟩.
Each many-body state |ψ⟩ = |config_up⟩ ⊗ |config_down⟩.
"""
function hamiltonian_element(basis::ManyBodyBasis,
                             idx_i::Int, idx_j::Int,
                             t::Float64, U::Float64)
    dim_up = length(basis.states_up)

    # Decompose indices into (up, down) configurations
    idx_i_up = ((idx_i - 1) % dim_up) + 1
    idx_i_down = ((idx_i - 1) ÷ dim_up) + 1
    idx_j_up = ((idx_j - 1) % dim_up) + 1
    idx_j_down = ((idx_j - 1) ÷ dim_up) + 1

    state_i_up = basis.states_up[idx_i_up]
    state_i_down = basis.states_down[idx_i_down]
    state_j_up = basis.states_up[idx_j_up]
    state_j_down = basis.states_down[idx_j_down]

    matrix_element = 0.0

    # Diagonal: interaction term U Σ_i n_{i,↑} n_{i,↓}
    if idx_i == idx_j
        for site in 1:basis.n_sites
            occ_up = is_occupied(state_i_up, site)
            occ_down = is_occupied(state_i_down, site)
            matrix_element += U * occ_up * occ_down
        end
    end

    # Off-diagonal: hopping terms
    # Spin-up hopping (down sector unchanged)
    if state_i_down == state_j_down
        for site in 1:basis.n_sites
            neighbors = get_neighbors(site, Lx, Ly)
            for neighbor in neighbors
                result = apply_hopping(state_j_up, site, neighbor)
                if result !== nothing
                    final_state, sign = result
                    if final_state == state_i_up
                        matrix_element += -t * sign
                    end
                end
            end
        end
    end

    # Spin-down hopping (up sector unchanged)
    if state_i_up == state_j_up
        for site in 1:basis.n_sites
            neighbors = get_neighbors(site, Lx, Ly)
            for neighbor in neighbors
                result = apply_hopping(state_j_down, site, neighbor)
                if result !== nothing
                    final_state, sign = result
                    if final_state == state_i_down
                        matrix_element += -t * sign
                    end
                end
            end
        end
    end

    return matrix_element
end

# =============================================================================
# 3. HAMILTONIAN MATRIX WRAPPER
# =============================================================================

"""
Wrapper to provide AbstractMatrix interface for the Hamiltonian.
"""
struct HubbardHamiltonian <: AbstractMatrix{Float64}
    basis::ManyBodyBasis
    t::Float64
    U::Float64
    dim::Int
end

Base.size(H::HubbardHamiltonian) = (H.dim, H.dim)
Base.eltype(::Type{HubbardHamiltonian}) = Float64

function Base.getindex(H::HubbardHamiltonian, i::Int, j::Int)
    return hamiltonian_element(H.basis, i, j, H.t, H.U)
end

# =============================================================================
# 4. BUILD HMATRIX (OPTIONAL - for large systems)
# =============================================================================

"""
Build HMatrix representation of Hamiltonian (useful for large Hilbert spaces).
For small systems, dense diagonalization may be faster.
"""
function build_hmatrix(H_op::HubbardHamiltonian; use_hmatrix=false)
    dim = H_op.dim

    println("\n" * "-"^70)
    println("Building Hamiltonian matrix...")
    println("-"^70)

    if !use_hmatrix || dim < 1000
        println("  Using dense matrix (dim = $dim)")
        return nothing  # Will use matrix-free approach
    end

    println("  Using HMatrix compression...")

    # Create trivial spatial coordinates (just use index as coordinate)
    X = [SVector{3}(Float64(i), 0.0, 0.0) for i in 1:dim]

    # Build cluster tree
    splitter = HMatrices.CardinalitySplitter(; nmax=50)
    Xclt = HMatrices.ClusterTree(X, splitter)
    Yclt = Xclt

    # Compression
    comp = HMatrices.PartialACA(; atol=1e-6, rank=30)

    # Assemble
    H = assemble_hmatrix(H_op, Xclt, Yclt; comp=comp)

    @printf("  Compression ratio: %.2f\n", compression_ratio(H))

    return H
end

# =============================================================================
# 5. DIAGONALIZATION WITH ARNOLDI
# =============================================================================

"""
Build sparse Hamiltonian matrix (much faster than matrix-free for repeated mat-vec).
"""
function build_sparse_hamiltonian(basis::ManyBodyBasis, t::Float64, U::Float64)
    dim_up = length(basis.states_up)
    dim_down = length(basis.states_down)
    dim = dim_up * dim_down

    println("\n" * "-"^70)
    println("Building sparse Hamiltonian matrix...")
    println("-"^70)
    @printf("  Hilbert space dimension: %d\n", dim)

    I_vals = Int[]
    J_vals = Int[]
    V_vals = Float64[]

    # Build sparse matrix elements
    @printf("  Computing matrix elements...")
    flush(stdout)

    for idx_i in 1:dim
        # Diagonal element
        idx_i_up = ((idx_i - 1) % dim_up) + 1
        idx_i_down = ((idx_i - 1) ÷ dim_up) + 1
        state_i_up = basis.states_up[idx_i_up]
        state_i_down = basis.states_down[idx_i_down]

        # Interaction term
        U_term = 0.0
        for site in 1:basis.n_sites
            occ_up = is_occupied(state_i_up, site)
            occ_down = is_occupied(state_i_down, site)
            U_term += U * occ_up * occ_down
        end
        if abs(U_term) > 1e-12
            push!(I_vals, idx_i)
            push!(J_vals, idx_i)
            push!(V_vals, U_term)
        end

        # Off-diagonal: hopping
        # Spin-up hopping
        for site in 1:basis.n_sites
            neighbors = get_neighbors(site, Lx, Ly)
            for neighbor in neighbors
                result = apply_hopping(state_i_up, site, neighbor)
                if result !== nothing
                    final_state, sign = result
                    # Find which basis state this corresponds to
                    idx_j_up = findfirst(==(final_state), basis.states_up)
                    if idx_j_up !== nothing
                        idx_j = idx_j_up + (idx_i_down - 1) * dim_up
                        push!(I_vals, idx_i)
                        push!(J_vals, idx_j)
                        push!(V_vals, -t * sign)
                    end
                end
            end
        end

        # Spin-down hopping
        for site in 1:basis.n_sites
            neighbors = get_neighbors(site, Lx, Ly)
            for neighbor in neighbors
                result = apply_hopping(state_i_down, site, neighbor)
                if result !== nothing
                    final_state, sign = result
                    idx_j_down = findfirst(==(final_state), basis.states_down)
                    if idx_j_down !== nothing
                        idx_j = idx_i_up + (idx_j_down - 1) * dim_up
                        push!(I_vals, idx_i)
                        push!(J_vals, idx_j)
                        push!(V_vals, -t * sign)
                    end
                end
            end
        end
    end

    @printf(" done\n")
    H_sparse = sparse(I_vals, J_vals, V_vals, dim, dim)
    @printf("  Sparse matrix: %d × %d\n", size(H_sparse, 1), size(H_sparse, 2))
    @printf("  Number of nonzero elements: %d\n", nnz(H_sparse))
    @printf("  Sparsity: %.2f%%\n", 100.0 * (1 - nnz(H_sparse) / dim^2))

    return H_sparse
end

"""
Find lowest eigenvalues using Arnoldi iteration.
"""
function diagonalize(H_sparse::SparseMatrixCSC)
    dim = size(H_sparse, 1)

    println("\n" * "-"^70)
    println("Diagonalizing with Arnoldi iteration...")
    println("-"^70)
    @printf("  Finding %d lowest eigenvalues\n", n_eigenvalues)

    # Matrix-vector product
    matvec = (v) -> H_sparse * v

    # Find smallest eigenvalues
    t_solve = @elapsed begin
        vals, vecs, info = eigsolve(matvec, dim, n_eigenvalues, :SR;
                                    issymmetric=true,
                                    krylovdim=max(30, 2*n_eigenvalues),
                                    maxiter=500,
                                    tol=1e-8)
    end

    @printf("  Converged: %d/%d eigenvalues\n", info.converged, n_eigenvalues)
    @printf("  Iterations: %d\n", info.numiter)
    @printf("  Time: %.3f seconds\n", t_solve)

    return real.(vals), vecs, info
end

# =============================================================================
# 6. ANALYSIS
# =============================================================================

"""
Analyze results.
"""
function analyze_results(eigenvalues, eigenvectors, basis)
    println("\n" * "="^70)
    println("RESULTS")
    println("="^70)

    println("\nLowest eigenvalues (energies in units of t):")
    for (i, E) in enumerate(eigenvalues[1:min(length(eigenvalues), n_eigenvalues)])
        @printf("  E_%d = %.8f eV = %.4f t\n", i-1, E, E/t)
    end

    E0 = eigenvalues[1]
    println("\nGround state energy:")
    @printf("  E₀ = %.8f eV\n", E0)
    @printf("  E₀/N = %.6f eV per site\n", E0 / N_sites)
    @printf("  E₀/t = %.6f\n", E0 / t)

    if length(eigenvalues) > 1
        gap = eigenvalues[2] - eigenvalues[1]
        @printf("\nEnergy gap to first excited state: Δ = %.6f eV\n", gap)
    end

    println("\nGround state wavefunction:")
    ψ0 = eigenvectors[1]
    @printf("  Norm: %.6f\n", norm(ψ0))
    @printf("  Largest amplitude: %.6f\n", maximum(abs.(ψ0)))

    # Find dominant configurations
    sorted_idx = sortperm(abs.(ψ0), rev=true)
    println("\n  Top 3 configurations:")
    for k in 1:min(3, length(ψ0))
        idx = sorted_idx[k]
        @printf("    |%d⟩: %.6f\n", idx-1, ψ0[idx])
    end

    return E0
end

# =============================================================================
# MAIN
# =============================================================================

function main()
    # Create basis
    basis = create_basis(N_sites, N_up, N_down)

    # Build sparse Hamiltonian
    H_sparse = build_sparse_hamiltonian(basis, t, U)

    # Diagonalize
    eigenvalues, eigenvectors, info = diagonalize(H_sparse)

    # Analyze
    E0 = analyze_results(eigenvalues, eigenvectors, basis)

    println("\n" * "="^70)
    println("Done!")
    println("="^70)

    return E0, eigenvalues, eigenvectors, basis
end

# Run
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
