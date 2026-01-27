#!/usr/bin/env julia
"""
Test get_reduced_grid function from Julia
"""

push!(LOAD_PATH, "/home/g/Research/FFirefly/src/module/imports")
include("/home/g/Research/FFirefly/src/module/imports/cpp_imports.jl")
using .Imports: get_reduced_grid
using Printf

function test_func(vec::Vector{Float64})
    """Test function that should be equal for symmetric points"""
    val = 0.0
    for v in vec
        val += cos(v)
    end
    return val
end

function ind_to_vec(inds::Vector{Int}, grid::Vector{Int})
    """Convert indices to momentum coordinates [-0.5, 0.5]"""
    v = Float64[]
    for (i, idx) in enumerate(inds)
        push!(v, -0.5 + idx / (grid[i] - 1))
    end
    return v
end

function test_symmetry_mapping(grid::Vector{Int}, lattice::String="SC")
    """Test that symmetry mapping produces valid equivalence classes"""
    println("\nTesting $(grid) grid with $(lattice) lattice:")

    reduced = get_reduced_grid(grid, lattice)

    println("  Number of symmetry groups: $(length(reduced))")

    # Count total points
    total_points = sum(length(group) for group in reduced)
    expected_points = prod(grid)

    println("  Total points: $(total_points) (expected: $(expected_points))")

    if total_points != expected_points
        println("  ❌ ERROR: Point count mismatch!")
        return false
    end

    # Test that points in each group have same test_func value
    all_valid = true
    for (i, group) in enumerate(reduced)
        values = Float64[]
        for point in group
            vec = ind_to_vec(point, grid)
            val = test_func(vec)
            push!(values, val)
        end

        # Check all values in group are equal (within tolerance)
        if length(values) > 1
            first_val = values[1]
            for val in values[2:end]
                if abs(val - first_val) > 1e-5
                    println("  ❌ Group $(i): Values not equal! $(values)")
                    all_valid = false
                    break
                end
            end
        end
    end

    if all_valid
        println("  ✓ All symmetry groups are valid!")
        return true
    else
        println("  ❌ Some groups have invalid symmetries")
        return false
    end
end

function print_grid_details(grid::Vector{Int}, lattice::String="SC")
    """Print detailed information about reduced grid"""
    println("\n" * "="^60)
    println("Detailed view of $(grid) grid with $(lattice) lattice:")
    println("="^60)

    reduced = get_reduced_grid(grid, lattice)

    for (i, group) in enumerate(reduced)
        println("\nGroup $(i-1): $(length(group)) points")
        for (j, point) in enumerate(group)
            vec = ind_to_vec(point, grid)
            val = test_func(vec)
            vec_str = join([Printf.@sprintf("%.3f", v) for v in vec], ", ")
            println("  $(point) -> [$(vec_str)] -> f=$(Printf.@sprintf("%.6f", val))")
        end
    end
end

function main()
    println("="^60)
    println("Testing get_reduced_grid from Julia")
    println("="^60)

    # Test different grid sizes
    test_cases = [
        ([5, 5], "SC"),
        ([6, 6], "SC"),
        ([4, 4], "SC"),
    ]

    all_passed = true
    for (grid, lattice) in test_cases
        passed = test_symmetry_mapping(grid, lattice)
        if !passed
            global all_passed = false
        end
    end

    # Print detailed view for one case
    print_grid_details([5, 5], "SC")

    println("\n" * "="^60)
    if all_passed
        println("✓ ALL TESTS PASSED")
    else
        println("❌ SOME TESTS FAILED")
    end
    println("="^60)
end

# Run tests
main()
