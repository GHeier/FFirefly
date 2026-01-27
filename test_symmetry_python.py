#!/usr/bin/env python3
"""
Test get_reduced_grid function from Python
"""
import sys
sys.path.insert(0, '/home/g/Research/FFirefly/src/module/imports')

from cpp_imports import get_reduced_grid
import math

def test_func(vec):
    """Test function that should be equal for symmetric points"""
    val = 0.0
    for v in vec:
        val += math.cos(v)
    return val

def ind_to_vec(inds, grid):
    """Convert indices to momentum coordinates [-0.5, 0.5]"""
    v = []
    for i, idx in enumerate(inds):
        v.append(-0.5 + idx / (grid[i] - 1))
    return v

def test_symmetry_mapping(grid, lattice="SC"):
    """Test that symmetry mapping produces valid equivalence classes"""
    print(f"\nTesting {grid} grid with {lattice} lattice:")

    reduced = get_reduced_grid(grid, lattice)

    print(f"  Number of symmetry groups: {len(reduced)}")

    # Count total points
    total_points = sum(len(group) for group in reduced)
    expected_points = 1
    for g in grid:
        expected_points *= g

    print(f"  Total points: {total_points} (expected: {expected_points})")

    if total_points != expected_points:
        print(f"  ❌ ERROR: Point count mismatch!")
        return False

    # Test that points in each group have same test_func value
    all_valid = True
    for i, group in enumerate(reduced):
        values = []
        for point in group:
            vec = ind_to_vec(point, grid)
            val = test_func(vec)
            values.append(val)

        # Check all values in group are equal (within tolerance)
        if len(values) > 1:
            first_val = values[0]
            for val in values[1:]:
                if abs(val - first_val) > 1e-5:
                    print(f"  ❌ Group {i}: Values not equal! {values}")
                    all_valid = False
                    break

    if all_valid:
        print(f"  ✓ All symmetry groups are valid!")
        return True
    else:
        print(f"  ❌ Some groups have invalid symmetries")
        return False

def print_grid_details(grid, lattice="SC"):
    """Print detailed information about reduced grid"""
    print(f"\n{'='*60}")
    print(f"Detailed view of {grid} grid with {lattice} lattice:")
    print(f"{'='*60}")

    reduced = get_reduced_grid(grid, lattice)

    for i, group in enumerate(reduced):
        print(f"\nGroup {i}: {len(group)} points")
        for j, point in enumerate(group):
            vec = ind_to_vec(point, grid)
            val = test_func(vec)
            print(f"  {point} -> {[f'{v:.3f}' for v in vec]} -> f={val:.6f}")

if __name__ == "__main__":
    print("="*60)
    print("Testing get_reduced_grid from Python")
    print("="*60)

    # Test different grid sizes
    test_cases = [
        ([5, 5], "SC"),
        ([6, 6], "SC"),
        ([4, 4], "SC"),
    ]

    all_passed = True
    for grid, lattice in test_cases:
        passed = test_symmetry_mapping(grid, lattice)
        if not passed:
            all_passed = False

    # Print detailed view for one case
    print_grid_details([5, 5], "SC")

    print("\n" + "="*60)
    if all_passed:
        print("✓ ALL TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
    print("="*60)
