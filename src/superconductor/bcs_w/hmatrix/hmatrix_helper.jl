module HMatrixHelper

export build_hmatrix

using HMatrices
using LinearAlgebra
using Printf
using StaticArrays

# Kernel matrix wrapper for k-k' points only (no frequency dependence)
struct KernelMatrixWrapper{F, T} <: AbstractMatrix{Float64}
    kernel::F
    kpoints::Vector{T}
    n_k::Int
end

function Base.getindex(K::KernelMatrixWrapper, i::Int, j::Int)
    k1 = K.kpoints[i]
    k2 = K.kpoints[j]
    return K.kernel(k1, k2, i, j)
end

Base.size(K::KernelMatrixWrapper) = (K.n_k, K.n_k)
Base.eltype(::Type{<:KernelMatrixWrapper}) = Float64

function build_spoints(kpoints)
    n_k = length(kpoints)
    k_dim = length(kpoints[1])

    # Build 3D clustering points from k-space
    if k_dim == 2
        X = [SVector{3}(Float64(kp[1]), Float64(kp[2]), 0.0) for kp in kpoints]
    elseif k_dim == 3
        X = [SVector{3}(Float64(kp[1]), Float64(kp[2]), Float64(kp[3])) for kp in kpoints]
    else
        error("K-points must be 2D or 3D, got dimension $k_dim")
    end
    return X
end

# Build HMatrix for k-k' kernel (no frequency dependence)
function build_hmatrix(kpoints, kernel_func; atol=1e-6, rank=20)
    n_k = length(kpoints)
    @printf("Building HMatrix for %d k-points\n", n_k)

    # Create kernel matrix wrapper
    K = KernelMatrixWrapper(kernel_func, kpoints, n_k)

    X = build_spoints(kpoints)

    # Build the hierarchical matrix
    println("Building HMatrix...")
    splitter = HMatrices.CardinalitySplitter(; nmax=50)
    Xclt = HMatrices.ClusterTree(X, splitter)
    Yclt = Xclt  # Same cluster tree for row/column spaces

    # Create compression method
    comp = HMatrices.PartialACA(; atol=atol, rank=rank)

    # Assemble HMatrix from matrix and cluster trees
    H = assemble_hmatrix(K, Xclt, Yclt; comp=comp)

    return H
end

end # module
