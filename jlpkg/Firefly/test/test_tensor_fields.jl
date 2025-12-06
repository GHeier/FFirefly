"""
Test script for tensor field save/load functionality in Julia.
Tests both 3D and 4D tensor fields.
"""

using Test
using Firefly

function test_tensor3_save_load()
    println("Testing 3D tensor save/load...")

    # Create test data: 3 spatial points, 2x2x2 tensor at each point
    num_tensors = 3
    ten_dim = 2
    mesh = Int32[3]
    domain = Float32[1.0;;]  # 1x1 matrix
    w_points = Float32[]

    # Create tensor data: shape (dim, dim, dim, num_tensors) in Julia column-major
    data = zeros(ComplexF32, ten_dim, ten_dim, ten_dim, num_tensors)
    for t in 1:num_tensors
        for i in 1:ten_dim
            for j in 1:ten_dim
                for k in 1:ten_dim
                    val = Float32(t + i + j + k - 4)  # Adjust for 1-based indexing
                    data[i, j, k, t] = ComplexF32(val, val / 10.0)
                end
            end
        end
    end

    # Save
    filename = "/tmp/test_tensor3_julia.h5"
    Firefly.Imports.save_data_tensor3(filename, data, num_tensors, ten_dim, true, mesh, domain, w_points)
    println("  Saved 3D tensor to $filename")

    # Load
    loaded = Firefly.Imports.BaseData(filename)
    println("  Loaded metadata: n_indices=$(loaded.n_indices), dim_indices=$(loaded.dim_indices)")

    # Check metadata
    @test loaded.n_indices == 3
    @test loaded.dim_indices == ten_dim
    @test loaded.nk == num_tensors

    # Get data back
    loaded_data = Firefly.Imports.get_data(loaded)
    println("  Loaded data shape: $(size(loaded_data))")

    # Check shape - Julia reshapes to (dim, dim, dim, num_tensors)
    @test size(loaded_data) == (ten_dim, ten_dim, ten_dim, num_tensors)

    # Check values
    max_diff = maximum(abs.(loaded_data .- data))
    println("  Max difference: $max_diff")
    @test max_diff < 1e-5

    # Cleanup
    rm(filename)
    println("  ✓ 3D tensor test passed!")
    return true
end

function test_tensor4_save_load()
    println("Testing 4D tensor save/load...")

    # Create test data: 4 spatial points, 2x2x2x2 tensor at each point
    num_tensors = 4
    ten_dim = 2
    mesh = Int32[2, 2]
    domain = Float32[1.0 0.0; 0.0 1.0]
    w_points = Float32[]

    # Create tensor data: shape (dim, dim, dim, dim, num_tensors) in Julia column-major
    data = zeros(ComplexF32, ten_dim, ten_dim, ten_dim, ten_dim, num_tensors)
    for t in 1:num_tensors
        for i in 1:ten_dim
            for j in 1:ten_dim
                for k in 1:ten_dim
                    for l in 1:ten_dim
                        val = Float32(t + i + j + k + l - 5)  # Adjust for 1-based indexing
                        data[i, j, k, l, t] = ComplexF32(val, val / 10.0)
                    end
                end
            end
        end
    end

    # Save
    filename = "/tmp/test_tensor4_julia.h5"
    Firefly.Imports.save_data_tensor4(filename, data, num_tensors, ten_dim, true, mesh, domain, w_points)
    println("  Saved 4D tensor to $filename")

    # Load
    loaded = Firefly.Imports.BaseData(filename)
    println("  Loaded metadata: n_indices=$(loaded.n_indices), dim_indices=$(loaded.dim_indices)")

    # Check metadata
    @test loaded.n_indices == 4
    @test loaded.dim_indices == ten_dim
    @test loaded.nk == num_tensors

    # Get data back
    loaded_data = Firefly.Imports.get_data(loaded)
    println("  Loaded data shape: $(size(loaded_data))")

    # Check shape
    @test size(loaded_data) == (ten_dim, ten_dim, ten_dim, ten_dim, num_tensors)

    # Check values
    max_diff = maximum(abs.(loaded_data .- data))
    println("  Max difference: $max_diff")
    @test max_diff < 1e-5

    # Cleanup
    rm(filename)
    println("  ✓ 4D tensor test passed!")
    return true
end

function test_tensor4_with_frequency()
    println("Testing 4D tensor with frequency dimension...")

    # Create test data: 2 w-points × 3 k-points, 2x2x2x2 tensor at each
    nw = 2
    nk = 3
    num_tensors = nw * nk
    ten_dim = 2
    mesh = Int32[3]
    domain = Float32[1.0;;]
    w_points = Float32[0.0, 1.0]

    # Create tensor data: shape (dim, dim, dim, dim, num_tensors)
    # In k-w ordering (C++ default): [k0w0, k0w1, k1w0, k1w1, k2w0, k2w1]
    data = zeros(ComplexF32, ten_dim, ten_dim, ten_dim, ten_dim, num_tensors)
    for k in 1:nk
        for w in 1:nw
            t = (k-1) * nw + w
            for i in 1:ten_dim
                for j in 1:ten_dim
                    for kk in 1:ten_dim
                        for l in 1:ten_dim
                            val = Float32(w + k + i + j + kk + l - 6)
                            data[i, j, kk, l, t] = ComplexF32(val, val / 10.0)
                        end
                    end
                end
            end
        end
    end

    # Save
    filename = "/tmp/test_tensor4_w_julia.h5"
    Firefly.Imports.save_data_tensor4(filename, data, num_tensors, ten_dim, true, mesh, domain, w_points)
    println("  Saved 4D tensor with frequency to $filename")

    # Load with k-w ordering (default)
    loaded = Firefly.Imports.BaseData(filename, "k-w")
    println("  Loaded metadata: n_indices=$(loaded.n_indices), nw=$(loaded.nw), nk=$(loaded.nk)")

    # Check metadata
    @test loaded.n_indices == 4
    @test loaded.nw == nw
    @test loaded.nk == nk

    # Get data back
    loaded_data = Firefly.Imports.get_data(loaded)
    println("  Loaded data shape: $(size(loaded_data))")

    # Check values
    max_diff = maximum(abs.(loaded_data .- data))
    println("  Max difference: $max_diff")
    @test max_diff < 1e-5

    # Cleanup
    rm(filename)
    println("  ✓ 4D tensor with frequency test passed!")
    return true
end

println("=" ^ 60)
println("Tensor Field Julia Interface Tests")
println("=" ^ 60)

@testset "Tensor Field Tests" begin
    test_tensor3_save_load()
    println()
    test_tensor4_save_load()
    println()
    test_tensor4_with_frequency()
end

println()
println("=" ^ 60)
println("✓ All Julia tensor field tests passed!")
println("=" ^ 60)
