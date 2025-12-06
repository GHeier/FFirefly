#!/usr/bin/env julia
"""
Example: BaseData and Field Usage in Julia

This script demonstrates:
1. Creating BaseData objects with various tensor ranks
2. Saving BaseData to HDF5
3. Loading BaseData from HDF5
4. Creating Field objects from data
5. Saving Field objects
6. Loading Field objects
7. Evaluating fields at arbitrary points
"""

# Add the Firefly package
push!(LOAD_PATH, joinpath(@__DIR__, "..", "jlpkg"))

using Firefly

function example_scalar_field()
    """Example 1: Scalar complex field with frequency dependence"""
    println("\n" * "="^60)
    println("Example 1: Scalar Complex Field (inds=[])")
    println("="^60)

    # Parameters
    nk = 20  # Number of k-points
    nw = 10  # Number of frequency points

    # Create scalar data (w-k ordering: frequency varies slowest)
    # Each point has a complex value
    data = ComplexF32[]
    for iw in 1:nw
        for ik in 1:nk
            # Value depends on position
            real_part = Float32((ik-1) / nk)
            imag_part = Float32((iw-1) / nw)
            push!(data, complex(real_part, imag_part))
        end
    end

    # Save using Firefly's save_data_scalar
    println("Creating scalar field...")
    mesh = Int32[nk]
    domain = Float32[1.0]  # 1D BZ spanning [-0.5, 0.5]
    w_points = Float32[(i-1) / nw for i in 1:nw]

    filename = "/tmp/example_scalar_field_jl.h5"
    save_data_scalar(filename, data, mesh, domain, w_points, is_complex=true)
    println("Saved BaseData to: $filename")

    # Load BaseData
    println("Loading BaseData...")
    basedata = BaseData(filename)
    println("  is_complex: $(basedata.is_complex)")
    println("  inds: $(basedata.inds)")
    println("  mesh: $(basedata.mesh)")
    println("  nk: $(basedata.nk)")
    println("  nw: $(basedata.nw)")

    # Create Field from file
    println("Creating Field_C from file...")
    field = Field_C(filename)

    # Evaluate at a point
    println("Evaluating field at k=0.0, w=0.5...")
    k = Vec(0.0)  # Center of BZ
    w = 0.5
    value = field(k, w)
    println("  Result: $value")

    # Save Field (should be identical to BaseData save)
    field_filename = "/tmp/example_scalar_field_from_field_jl.h5"
    save(field, field_filename)
    println("Saved Field to: $field_filename")

    # Load Field again
    field2 = Field_C(field_filename)
    value2 = field2(k, w)
    println("Loaded Field evaluation: $value2")
    println("Values match: $(abs(value - value2) < 1e-6)")
end


function example_matrix_field()
    """Example 2: Matrix field (rank-2 tensor)"""
    println("\n" * "="^60)
    println("Example 2: Complex Matrix Field (inds=[3, 3])")
    println("="^60)

    # Parameters
    nk = 10  # Total k-points (10x10 mesh)
    mat_dim = 3  # 3x3 matrices

    # Create matrix data: vector of 3x3 matrices
    # Each matrix at each k-point
    data = Vector{Matrix{ComplexF32}}(undef, nk * nk)
    for ik in 1:(nk * nk)
        matrix = Matrix{ComplexF32}(undef, mat_dim, mat_dim)
        for i in 1:mat_dim
            for j in 1:mat_dim
                # Matrix elements depend on position and indices
                val = Float32((ik-1) / (nk*nk) + (i + j - 2) / 10.0)
                matrix[i, j] = complex(val, val / 10.0)
            end
        end
        data[ik] = matrix
    end

    # Save matrix field
    filename = "/tmp/example_matrix_field_jl.h5"
    mesh = Int32[nk, nk]
    domain = Float32[1.0 0.0; 0.0 1.0]  # 2D BZ
    w_points = Float32[]
    inds = Int32[mat_dim, mat_dim]

    println("Creating $(mat_dim)x$(mat_dim) matrix field...")
    save_data_matrix(filename, data, mesh, domain, w_points, inds, is_complex=true)
    println("Saved to: $filename")

    # Load and inspect BaseData
    println("Loading BaseData...")
    basedata = BaseData(filename)
    println("  inds: $(basedata.inds)")
    println("  rank: $(length(basedata.inds))")
    println("  dimension: $(basedata.dimension)")

    # Create Field_CM
    println("Creating Field_CM...")
    field = Field_CM(filename)

    # Evaluate at a point
    k = Vec(0.1, 0.2)
    println("Evaluating at k=($(k.x), $(k.y))...")
    matrix = field(k)
    println("  Result shape: $(size(matrix))")
    println("  Matrix[1,1] = $(matrix[1,1])")
    println("  Matrix[2,3] = $(matrix[2,3])")

    # Save and reload Field
    field_file = "/tmp/example_matrix_field_from_field_jl.h5"
    save(field, field_file)
    field2 = Field_CM(field_file)
    matrix2 = field2(k)
    println("Round-trip match: $(abs(matrix[1,1] - matrix2[1,1]) < 1e-6)")
end


function example_nonuniform_matrix()
    """Example 3: Non-uniform matrix (2x3)"""
    println("\n" * "="^60)
    println("Example 3: Non-Uniform Matrix Field (inds=[2, 3])")
    println("="^60)

    # Create 2x3 matrix data
    nk = 15
    dim1, dim2 = 2, 3

    data = Vector{Matrix{ComplexF32}}(undef, nk)
    for ik in 1:nk
        matrix = Matrix{ComplexF32}(undef, dim1, dim2)
        for i in 1:dim1
            for j in 1:dim2
                val = Float32(((ik-1) * dim1 * dim2 + (i-1) * dim2 + (j-1)) / (nk * dim1 * dim2))
                matrix[i, j] = complex(val, -val)
            end
        end
        data[ik] = matrix
    end

    filename = "/tmp/example_nonuniform_matrix_jl.h5"
    mesh = Int32[nk]
    domain = Float32[1.0]
    w_points = Float32[]
    inds = Int32[dim1, dim2]

    println("Creating $(dim1)x$(dim2) non-uniform matrix field...")
    save_data_matrix(filename, data, mesh, domain, w_points, inds, is_complex=true)

    # Load and use
    field = Field_CM(filename)
    k = Vec(0.0)
    result = field(k)
    println("  Result shape: $(size(result))")
    println("  Non-square matrix verified!")
end


function example_4d_vertex()
    """Example 4: 4D vertex tensor (multi-orbital interaction)"""
    println("\n" * "="^60)
    println("Example 4: 4D Vertex Tensor (inds=[2, 2, 2, 2])")
    println("="^60)

    # Create 2x2x2x2 tensor data (2-orbital vertex)
    nk = 25  # 5x5 k-mesh
    dim = 2
    nw = 3

    # Create 4D tensor data
    data = Vector{Array{ComplexF32, 4}}(undef, nw * nk)
    idx_global = 1
    for iw in 1:nw
        for ik in 1:nk
            tensor = Array{ComplexF32, 4}(undef, dim, dim, dim, dim)
            for i in 1:dim
                for j in 1:dim
                    for k in 1:dim
                        for l in 1:dim
                            idx = (i-1)*dim^3 + (j-1)*dim^2 + (k-1)*dim + (l-1)
                            val = Float32((ik + iw - 2) / (nk + nw) + idx / 100.0)
                            tensor[i, j, k, l] = complex(val, val/20.0)
                        end
                    end
                end
            end
            data[idx_global] = tensor
            idx_global += 1
        end
    end

    filename = "/tmp/example_4d_vertex_jl.h5"
    mesh = Int32[5, 5]
    domain = Float32[6.28 0.0; 0.0 6.28]  # Full BZ
    w_points = Float32[-0.1, 0.0, 0.1]
    inds = Int32[dim, dim, dim, dim]

    println("Creating $(dim)x$(dim)x$(dim)x$(dim) vertex tensor...")
    save_data_tensor4(filename, data, mesh, domain, w_points, inds, is_complex=true)

    # Load as Field_CM
    field = Field_CM(filename)
    k = Vec(3.14, 3.14)  # BZ center
    w = 0.0

    # Field_CM flattens 4D tensor to [d1*d2][d3*d4] matrix
    result = field(k, w)
    println("  Result shape (flattened): $(size(result))")
    println("  Expected: [$(dim*dim)][$(dim*dim)] = [4][4]")
    println("  result[1,1] = $(result[1,1])")
end


function example_vector_field()
    """Example 5: Vector field (rank-1)"""
    println("\n" * "="^60)
    println("Example 5: Vector Field (inds=[4])")
    println("="^60)

    # Create 4-component vector field
    nk = 10
    vec_dim = 4

    data = Vector{Vector{ComplexF32}}(undef, nk)
    for ik in 1:nk
        vector = Vector{ComplexF32}(undef, vec_dim)
        for i in 1:vec_dim
            val = Float32((ik-1) / nk + (i-1) / vec_dim)
            vector[i] = complex(val, -val/2.0)
        end
        data[ik] = vector
    end

    filename = "/tmp/example_vector_field_jl.h5"
    mesh = Int32[nk]
    domain = Float32[1.0]
    w_points = Float32[]
    inds = Int32[vec_dim]

    println("Creating $(vec_dim)-component vector field...")
    save_data_vector(filename, data, mesh, domain, w_points, inds, is_complex=true)

    # Load as Field_CM
    # For rank-1, Field_CM wraps as column matrix [vec_dim][1]
    field = Field_CM(filename)
    k = Vec(0.5)
    result = field(k)

    println("  Result shape (as column matrix): $(size(result))")
    println("  Vector components:")
    for i in 1:size(result, 1)
        println("    [$i]: $(result[i,1])")
    end
end


function main()
    """Run all examples"""
    println("\n" * "="^60)
    println("BaseData and Field Examples - Julia")
    println("="^60)

    example_scalar_field()
    example_matrix_field()
    example_nonuniform_matrix()
    example_4d_vertex()
    example_vector_field()

    println("\n" * "="^60)
    println("All examples completed successfully!")
    println("="^60)
    println("\nFiles created in /tmp/:")
    println("  - example_scalar_field_jl.h5")
    println("  - example_matrix_field_jl.h5")
    println("  - example_nonuniform_matrix_jl.h5")
    println("  - example_4d_vertex_jl.h5")
    println("  - example_vector_field_jl.h5")
    println("\nYou can inspect these files with h5dump or h5ls")
end


# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
