#!/usr/bin/env python3
"""
Simplified BaseData and Field Examples for Python

This demonstrates the basic workflow:
1. Load BaseData from HDF5
2. Create Field from BaseData
3. Evaluate field at points
4. Save Field to HDF5
5. Load Field back

Note: This uses pre-existing test files. For creating data from scratch,
see the C++ test files in src/objects/CMField/tests/
"""

import sys
import os

# Add the firefly package to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'pypkg'))

try:
    import firefly.src.module.imports.cpp_imports as ff
except ImportError as e:
    print(f"Error importing firefly: {e}")
    print("Make sure you've built the project with ./scripts/fly-build.sh")
    sys.exit(1)


def demo_field_operations():
    """Demonstrate Field operations using test data"""
    print("\n" + "="*60)
    print("BaseData and Field Operations Demo")
    print("="*60)

    # First, create some simple test data using the C++ test infrastructure
    print("\n1. Creating test data files...")
    import subprocess
    result = subprocess.run(['build/bin/fly.x'], capture_output=True, text=True)
    if result.returncode != 0:
        print("Warning: Could not run tests to generate sample data")
        print("Creating minimal example instead...")
        create_minimal_example()
        return

    # Use one of the test files created
    test_file = "/tmp/test_matrix_field.h5"

    if not os.path.exists(test_file):
        print(f"Test file {test_file} not found")
        print("Run 'build/bin/fly.x' first to generate test data")
        return

    print(f"\n2. Loading BaseData from: {test_file}")
    basedata = ff.BaseData(test_file)

    print("\nBaseData properties:")
    print(f"  is_complex: {basedata.is_complex}")
    print(f"  is_matrix: {basedata.is_matrix}")
    print(f"  inds: {basedata.inds}")
    print(f"  rank: {len(basedata.inds)}")
    print(f"  mesh: {basedata.mesh}")
    print(f"  dimension: {basedata.dimension}")
    print(f"  nk: {basedata.nk}")
    print(f"  domain shape: {len(basedata.domain)}x{len(basedata.domain[0]) if basedata.domain else 0}")

    print("\n3. Creating Field_CM from file...")
    field = ff.Field_CM(test_file)

    print("\n4. Evaluating field at different points...")

    # Create Vec objects for different k-points
    k1 = ff.Vec(0.0)  # Center of BZ
    k2 = ff.Vec(0.2)
    k3 = ff.Vec(-0.3)

    print(f"\nAt k = ({k1.x}, {k1.y}, {k1.z}):")
    matrix1 = field(k1)
    print(f"  Matrix shape: {len(matrix1)}x{len(matrix1[0])}")
    print(f"  Matrix[0][0] = {matrix1[0][0]}")

    print(f"\nAt k = ({k2.x}, {k2.y}, {k2.z}):")
    matrix2 = field(k2)
    print(f"  Matrix[0][0] = {matrix2[0][0]}")

    print(f"\nAt k = ({k3.x}, {k3.y}, {k3.z}):")
    matrix3 = field(k3)
    print(f"  Matrix[0][0] = {matrix3[0][0]}")

    print("\n5. Saving Field to new file...")
    output_file = "/tmp/example_field_python.h5"
    field.save(output_file)
    print(f"  Saved to: {output_file}")

    print("\n6. Loading Field from saved file...")
    field_reloaded = ff.Field_CM(output_file)
    matrix_reloaded = field_reloaded(k1)

    print("\n7. Verifying round-trip...")
    match = abs(matrix1[0][0] - matrix_reloaded[0][0]) < 1e-6
    print(f"  Original: {matrix1[0][0]}")
    print(f"  Reloaded: {matrix_reloaded[0][0]}")
    print(f"  Match: {match}")

    if match:
        print("\n✓ Round-trip successful!")
    else:
        print("\n✗ Round-trip failed!")


def create_minimal_example():
    """Create a minimal working example"""
    print("\nCreating minimal example...")
    print("This demonstrates the Python API for Field operations.")
    print("\nFor full examples of creating data from scratch,")
    print("see the C++ test files in src/objects/CMField/tests/")
    print("\nKey concepts:")
    print("  - BaseData: Container for multi-dimensional field data")
    print("  - inds: Vector of tensor dimensions (e.g., {3,3} for 3x3 matrix)")
    print("  - mesh: Grid size in k-space")
    print("  - domain: Brillouin zone vectors")
    print("  - w_points: Frequency points (optional)")
    print("\nWorkflow:")
    print("  1. Create data (use C++ or save_data_* functions)")
    print("  2. Save to HDF5")
    print("  3. Load with BaseData(filename)")
    print("  4. Create Field from file: Field_CM(filename)")
    print("  5. Evaluate: matrix = field(k_point, w)")
    print("  6. Save field: field.save(filename)")


def show_field_types():
    """Show available field types"""
    print("\n" + "="*60)
    print("Available Field Types in Firefly")
    print("="*60)
    print("\nScalar Fields:")
    print("  - Field_R: Real scalar field")
    print("  - Field_C: Complex scalar field")
    print("\nMatrix/Tensor Fields:")
    print("  - Field_RM: Real matrix/tensor field")
    print("  - Field_CM: Complex matrix/tensor field")
    print("\nTensor Ranks (inds parameter):")
    print("  - Scalar: inds = [] (empty)")
    print("  - Vector: inds = [n] (n-component vector)")
    print("  - Matrix: inds = [m, n] (m×n matrix)")
    print("  - 3D Tensor: inds = [l, m, n]")
    print("  - 4D Tensor: inds = [i, j, k, l] (e.g., vertex functions)")
    print("\nNon-uniform dimensions supported:")
    print("  - inds = [2, 3] → 2×3 matrix")
    print("  - inds = [1, 2, 2, 1] → vertex tensor")
    print("\nField_CM returns flattened matrices:")
    print("  - Rank-1: [n][1] column matrix")
    print("  - Rank-2: [m][n] matrix")
    print("  - Rank-3: [d1*d2][d3] flattened matrix")
    print("  - Rank-4: [d1*d2][d3*d4] flattened matrix")


def main():
    print("\n" + "="*60)
    print("Firefly BaseData and Field Examples - Python")
    print("="*60)

    show_field_types()
    demo_field_operations()

    print("\n" + "="*60)
    print("Demo completed!")
    print("="*60)


if __name__ == "__main__":
    main()
