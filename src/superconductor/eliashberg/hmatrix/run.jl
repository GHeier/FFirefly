using Firefly
cfg = Firefly.Config

using HMatrices
using KrylovKit
using LinearAlgebra
using Printf
using StaticArrays
# Load relevant variables from the configuration
kmesh = cfg.k_mesh


function load_surface()
    eps_func = (k) -> Firefly.epsilon(1, [Float64(k.x), Float64(k.y), Float64(k.z)])

    # Create Surface at the Fermi energy (mu)
    surf = Firefly.Surface(eps_func, Float32(mu))
    faces = get_faces(surf)
    println("Found $(length(faces)) k-points on Fermi surface")

    return faces
end

function kernel(k1, k2)
    dk = k1 .- k2
    q = LinearAlgebra.norm(dk)
    # Regularized Coulomb: 1/(q^2 + λ^2)
    λ = 0.1
    return 1.0 / (q^2 + λ^2)
end

struct KernelMatrixWrapper{F, T} <: AbstractMatrix{Float64}
    kernel::F
    kpoints::Vector{T}
end

Base.getindex(K::KernelMatrixWrapper, i::Int, j::Int) = K.kernel(K.kpoints[i], K.kpoints[j])
Base.size(K::KernelMatrixWrapper) = (length(K.kpoints), length(K.kpoints))
Base.eltype(::Type{<:KernelMatrixWrapper}) = Float64

# Build HMatrix
function build_hmatrix(kpoints, kernel_func; atol=1e-6, rank=20)
    n = length(kpoints)

    # Convert kpoints to vector of SVectors for HMatrices
    X = [SVector{3}(kp) for kp in kpoints]

    # Create an indexable kernel matrix wrapper
    K = KernelMatrixWrapper(kernel_func, kpoints)

    # Build the hierarchical matrix
    println("Building HMatrix...")
    splitter = HMatrices.CardinalitySplitter(; nmax=50)
    Xclt = HMatrices.ClusterTree(X, splitter)
    Yclt = Xclt  # Same cluster tree for symmetric matrix

    # Create compression method
    comp = HMatrices.PartialACA(; atol=atol, rank=rank)

    # Assemble HMatrix from matrix and cluster trees
    H = assemble_hmatrix(K, Xclt, Yclt; comp=comp)

    return H
end

function run():
    # Main function call goes here
    println("="^60)
    println("HMatrix + Lanczos Eigenvalue Solver Test")
    println("="^60)

    # Parameters
    mu = -1.0  # Fermi energy
    println("\nGenerating k-points from Fermi surface at μ = $mu")
    println("Testing epsilon at a few k-points:")
    for ktest in [[0.0, 0.0, 0.0], [π, π, π], [0.5, 0.5, 0.5]]
        ε = Firefly.epsilon(1, ktest)
        println("  ε($(ktest)) = $ε")
    end

    kpoints = generate_kpoints_from_surface(mu)
    n = length(kpoints)
    println("Total k-points: $n")

    # Build HMatrix
    println("\n" * "="^60)
    t_hmat = @elapsed H = build_hmatrix(kpoints, kernel; atol=1e-4, rank=30)
    @printf("HMatrix built in %.3f seconds\n", t_hmat)
    @printf("Compression ratio: %.2f\n", compression_ratio(H))
    @printf("Memory savings: %.2f%%\n", (1 - compression_ratio(H)) * 100)

    # Solve with HMatrix using Lanczos
    println("\n" * "="^60)
    println("Solving HMatrix with Lanczos (KrylovKit)...")

    # Create a function that applies H to a vector using mul!
    # HMatrices.jl supports matrix-vector multiplication via mul!
    hmv = (v) -> begin
        result = similar(v)
        mul!(result, H, v)
        return result
    end

    # Use the matrix-vector product function with eigsolve
    t_lanczos_hmat = @elapsed begin
        vals_hmat, vecs_hmat, info_hmat = eigsolve(hmv, n, 1, :LM;
                                                    issymmetric=true,
                                                    krylovdim=30,
                                                    maxiter=200,
                                                    tol=1e-8)
    end
    @printf("Time: %.3f seconds\n", t_lanczos_hmat)
    @printf("Largest eigenvalue: %.10f\n", real(vals_hmat[1]))
    @printf("Converged: %s (iterations: %d)\n", info_hmat.converged > 0, info_hmat.numiter)

    return real(vals_hmat[1])
end

if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



