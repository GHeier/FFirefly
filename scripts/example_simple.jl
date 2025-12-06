#!/usr/bin/env julia
"""
BaseData and Fields Tutorial. Very simple, just go to bottom of file
"""

using LinearAlgebra
using Printf

using Firefly
using Firefly.Imports

const file = "temp_data.h5"

# f(w,k) = exp(-k^2/(0.05+w)^2)
function gaussian(k::Vector{<:Real}, w::Real)
    kmag = k[1]^2 + k[2]^2
    sigma = 0.2 + w
    return exp(-kmag / sigma^2)
end

function get_data()
    # w-points are defined, for this example k-points are assumed to lie on a mesh
    # If non-mesh functions are desired, they can also be passed in as a list of points
    w_points = [0.0, 0.5, 1.0]
    nw = length(w_points)
    nk = 100
    k_mesh = [nk, nk]  # 2D Mesh
    println("w-points: ", w_points)
    println("k-mesh: ", k_mesh)

    # K-space mesh
    kx = range(-0.5, 0.5, length=nk)
    ky = range(-0.5, 0.5, length=nk)

    # Create k-point grid
    kgrid = [[x, y] for x in kx, y in ky]
    kgrid_flat = vec(kgrid)  # Flatten to 1D array of vectors

    # Brillouin Zone (in 3D by default)
    BZ = Float32[2*π 0.0 0.0; 0.0 2*π 0.0; 0.0 0.0 2*π]
    # Brillouin Zone changed to 2D
    BZ = BZ[1:2, 1:2]

    # Transform k-points
    kpts = [vec(Float32.(BZ * k)) for k in kgrid_flat]
    kpts_mat = reduce(hcat, kpts)'  # Convert to matrix for min/max
    println("kmin: ", minimum(kpts_mat, dims=1))
    println("kmax: ", maximum(kpts_mat, dims=1))
    println("K-point list shape: ", size(kpts_mat))

    # Create Gaussian data: f(w,k) = exp(-k^2/(0.05+w)^2)
    # w-k ordering: frequency varies slowest
    println("\nCreating Gaussian data on $(nw)x$(nk) grid...")

    # Python approach: data[i, :, :] = gaussian(kpts.T, w).reshape(nk, nk) with C-order
    # This fills the array row-by-row in row-major order
    # Julia's reshape fills column-by-column in column-major order
    # So Julia's reshape creates the transpose of what Python creates
    # Solution: Don't reshape - build the array directly in the correct order
    data = zeros(Float32, nw, nk, nk)
    for i in 1:nw
        # Build data in row-major order to match Python
        for ix in 1:nk
            for iy in 1:nk
                idx = (ix-1)*nk + iy  # Row-major index
                data[i, ix, iy] = Float32(gaussian(kpts[idx], w_points[i]))
            end
        end
    end
    println("data shape: ", size(data))
    return data, k_mesh, BZ, w_points
end

function BaseData_example(data, k_mesh, BZ, w_points)
    # Step 1: Save using save_data!
    println("\nStep 1: Save data")
    save_data!(file, data, mesh=k_mesh, domain=BZ, w_points=w_points)
    println("  Saved to: ", file)

    # Step 2: Load into BaseData
    println("\nStep 2: Load into BaseData")
    bd = Imports.BaseData(file)
    println("  Loaded BaseData:")
    println("    is_complex: ", bd.is_complex)
    println("    inds: ", bd.inds)
    println("    rank: ", length(bd.inds))
    println("    mesh: ", bd.mesh)
    println("    nk: ", bd.nk)
    println("    nw: ", bd.nw)

    # Access data through the data property
    println("\n  Accessing data through bd data:")
    bd_data = Imports.get_data(bd)
    println("    data shape: ", size(bd_data))
    println("    data type: ", eltype(bd_data))
    @printf("    data min/max: %.6f / %.6f\n", minimum(bd_data), maximum(bd_data))

    # Step 3: Save BaseData
    println("\nStep 3: Save BaseData")
    Imports.save!(bd, file)
    return bd
end

function Field_R_example()
    println("\nLoad into Field_R (Real Scalar Field) from file")
    field = Field_R(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print("    Evaluating at k=(0.0,0.0), w=0.5:")
    @printf("    Expected, found: %.6f, %f\n", expected, val)

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5]
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print("    Evaluating at k=(-0.2, 0.5), w=1.0:")
    println("    Expected, found: ", expected, ", ", val)
end

function Field_C_example()
    println("\nLoad into Field_C (Complex Scalar Field) from file")
    field = Field_C(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print("    Evaluating at k=(0.0,0.0), w=0.5:")
    println("    Expected, found: ", expected, ", ", val)

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5]
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print("    Evaluating at k=(-0.2, 0.5), w=1.0:")
    println("    Expected, found: ", expected, ", ", val)
end

function Field_CM_example()
    println("\nLoad into Field_CM (Complex Matrix Field) from file")
    field = Field_CM(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print("    Evaluating at k=(0.0,0.0), w=0.5:")
    println("    Expected, found: ", expected, ", ", val)

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5]
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print("    Evaluating at k=(-0.2, 0.5), w=1.0:")
    println("    Expected, found: ", expected, ", ", val)
end

function Field_RM_example()
    println("\nLoad into Field_RM (Real Matrix Field) from file")
    field = Field_RM(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print("    Evaluating at k=(0.0,0.0), w=0.5:")
    println("    Expected, found: ", expected, ", ", val)

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5]
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print("    Evaluating at k=(-0.2, 0.5), w=1.0:")
    println("    Expected, found: ", expected, ", ", val)
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    # Create data on a 3D w-k grid
    data, k_mesh, BZ, w_points = get_data()
    save_data!(file, data, mesh=k_mesh, domain=BZ, w_points=w_points)

    bd = BaseData_example(data, k_mesh, BZ, w_points)
    # Note: file already saved in BaseData_example, no need to save again

    # Make a Real Scalar Field
    Field_R_example()

    # Make Data Complex
    bd_data = Imports.get_data(bd)
    bd._data = Complex{Float32}.(bd_data)
    Imports.save!(bd, file)
    # Make a Complex Scalar Field
    Field_C_example()

    # Reshape to matrix format
    bd_data = Imports.get_data(bd)
    # get_data now returns data in column-major order, ready for Julia reshape
    bd._data = reshape(bd_data, bd.nw, bd.mesh[1], bd.mesh[2], 1, 1)
    bd.inds = [1, 1]
    Imports.save!(bd, file)

    Field_CM_example()

    # Convert to real
    bd_data = Imports.get_data(bd)
    # bd_data is already in correct Julia shape from bd._data, just convert type
    bd._data = Float32.(real.(bd_data))
    Imports.save!(bd, file)

    Field_RM_example()

    if isfile(file)
        rm(file)
        println("\nFile '", file, "' deleted successfully.")
    end
end
