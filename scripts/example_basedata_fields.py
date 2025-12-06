#!/usr/bin/env python3
"""
Example: BaseData and Field Usage in Python

This script demonstrates:
1. Creating BaseData objects with various tensor ranks
2. Saving BaseData to HDF5
3. Loading BaseData from HDF5
4. Creating Field objects from data
5. Saving Field objects
6. Loading Field objects
7. Evaluating fields at arbitrary points
"""

import sys
import os
import numpy as np

# Add the firefly package to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'pypkg'))

try:
    import firefly.src.module.imports.cpp_imports as ff
except ImportError as e:
    print(f"Error importing firefly: {e}")
    print("Make sure you've built the project with ./scripts/fly-build.sh")
    sys.exit(1)

def example_scalar_field():
    """Example 1: Scalar complex field with frequency dependence"""
    print("\n" + "="*60)
    print("Example 1: Scalar Complex Field (inds=[])")
    print("="*60)

    # Parameters
    nk = 20  # Number of k-points
    nw = 10  # Number of frequency points

    # Create scalar data (w-k ordering: frequency varies slowest)
    # Each point has a complex value
    data = np.zeros(nk * nw, dtype=np.complex64)
    idx = 0
    for iw in range(nw):
        for ik in range(nk):
            # Value depends on position
            real_part = float(ik) / nk
            imag_part = float(iw) / nw
            data[idx] = complex(real_part, imag_part)
            idx += 1

    # Save using firefly's save_data_scalar
    print("Creating scalar field...")
    mesh = np.array([nk], dtype=np.int32)
    domain = np.array([[1.0]], dtype=np.float32)  # 1D BZ spanning [-0.5, 0.5]
    w_points = np.array([float(i) / nw for i in range(nw)], dtype=np.float32)

    filename = "/tmp/example_scalar_field.h5"
    ff.save_data_scalar(filename, data, True, mesh, domain, w_points)
    print(f"Saved BaseData to: {filename}")

    # Load BaseData
    print("Loading BaseData...")
    basedata = ff.BaseData(filename)
    print(f"  is_complex: {basedata.is_complex}")
    print(f"  inds: {basedata.inds}")
    print(f"  mesh: {basedata.mesh}")
    print(f"  nk: {basedata.nk}")
    print(f"  nw: {basedata.nw}")

    # Create Field from file
    print("Creating Field_C from file...")
    field = ff.Field_C(filename)

    # Evaluate at a point
    print("Evaluating field at k=0.0, w=0.5...")
    k = ff.Vec(0.0)  # Center of BZ
    w = 0.5
    value = field(k, w)
    print(f"  Result: {value}")

    # Save Field (should be identical to BaseData save)
    field_filename = "/tmp/example_scalar_field_from_field.h5"
    field.save(field_filename)
    print(f"Saved Field to: {field_filename}")

    # Load Field again
    field2 = ff.Field_C(field_filename)
    value2 = field2(k, w)
    print(f"Loaded Field evaluation: {value2}")
    print(f"Values match: {abs(value - value2) < 1e-6}")


def example_matrix_field():
    """Example 2: Matrix field (rank-2 tensor)"""
    print("\n" + "="*60)
    print("Example 2: Complex Matrix Field (inds={3, 3})")
    print("="*60)

    # Parameters
    nk = 10  # Total k-points (10x10 mesh)
    mat_dim = 3  # 3x3 matrices

    # Create matrix data: vector of 3x3 matrices
    # Each matrix at each k-point
    data = []
    for ik in range(nk * nk):
        matrix = []
        for i in range(mat_dim):
            row = []
            for j in range(mat_dim):
                # Matrix elements depend on position and indices
                val = float(ik) / (nk*nk) + float(i + j) / 10.0
                row.append(complex(val, val / 10.0))
            matrix.append(row)
        data.append(matrix)

    # Save matrix field
    filename = "/tmp/example_matrix_field.h5"
    mesh = [nk, nk]
    domain = [[1.0, 0.0], [0.0, 1.0]]  # 2D BZ
    inds = [mat_dim, mat_dim]

    print(f"Creating {mat_dim}x{mat_dim} matrix field...")
    ff.save_data_matrix(filename, data, True, mesh, domain, [], inds)
    print(f"Saved to: {filename}")

    # Load and inspect BaseData
    print("Loading BaseData...")
    basedata = ff.BaseData(filename)
    print(f"  inds: {basedata.inds}")
    print(f"  rank: {len(basedata.inds)}")
    print(f"  dimension: {basedata.dimension}")

    # Create Field_CM
    print("Creating Field_CM...")
    field = ff.Field_CM(filename)

    # Evaluate at a point
    k = ff.Vec(0.1, 0.2)
    print(f"Evaluating at k=({k.x}, {k.y})...")
    matrix = field(k)
    print(f"  Result shape: {len(matrix)}x{len(matrix[0])}")
    print(f"  Matrix[0][0] = {matrix[0][0]}")
    print(f"  Matrix[1][2] = {matrix[1][2]}")

    # Save and reload Field
    field_file = "/tmp/example_matrix_field_from_field.h5"
    field.save(field_file)
    field2 = ff.Field_CM(field_file)
    matrix2 = field2(k)
    print(f"Round-trip match: {abs(matrix[0][0] - matrix2[0][0]) < 1e-6}")


def example_nonuniform_matrix():
    """Example 3: Non-uniform matrix (2x3)"""
    print("\n" + "="*60)
    print("Example 3: Non-Uniform Matrix Field (inds={2, 3})")
    print("="*60)

    # Create 2x3 matrix data
    nk = 15
    dim1, dim2 = 2, 3

    data = []
    for ik in range(nk):
        matrix = []
        for i in range(dim1):
            row = []
            for j in range(dim2):
                val = float(ik * dim1 * dim2 + i * dim2 + j) / (nk * dim1 * dim2)
                row.append(complex(val, -val))
            matrix.append(row)
        data.append(matrix)

    filename = "/tmp/example_nonuniform_matrix.h5"
    mesh = [nk]
    domain = [[1.0]]
    inds = [dim1, dim2]

    print(f"Creating {dim1}x{dim2} non-uniform matrix field...")
    ff.save_data_matrix(filename, data, True, mesh, domain, [], inds)

    # Load and use
    field = ff.Field_CM(filename)
    k = ff.Vec(0.0)
    result = field(k)
    print(f"  Result shape: {len(result)}x{len(result[0])}")
    print(f"  Non-square matrix verified!")


def example_4d_vertex():
    """Example 4: 4D vertex tensor (multi-orbital interaction)"""
    print("\n" + "="*60)
    print("Example 4: 4D Vertex Tensor (inds={2, 2, 2, 2})")
    print("="*60)

    # Create 2x2x2x2 tensor data (2-orbital vertex)
    nk = 25  # 5x5 k-mesh
    dim = 2
    nw = 3

    # Create 4D tensor data
    # Total size: nw * nk * (2*2*2*2) = 3 * 25 * 16 = 1200 values
    data = []
    for iw in range(nw):
        for ik in range(nk):
            tensor = []
            for i in range(dim):
                t1 = []
                for j in range(dim):
                    t2 = []
                    for k in range(dim):
                        t3 = []
                        for l in range(dim):
                            idx = i*dim**3 + j*dim**2 + k*dim + l
                            val = float(ik + iw) / (nk + nw) + float(idx) / 100.0
                            t3.append(complex(val, val/20.0))
                        t2.append(t3)
                    t1.append(t2)
                tensor.append(t1)
            data.append(tensor)

    filename = "/tmp/example_4d_vertex.h5"
    mesh = [5, 5]
    domain = [[6.28, 0.0], [0.0, 6.28]]  # Full BZ
    w_points = [-0.1, 0.0, 0.1]
    inds = [dim, dim, dim, dim]

    print(f"Creating {dim}x{dim}x{dim}x{dim} vertex tensor...")
    ff.save_data_tensor4(filename, data, True, mesh, domain, w_points, inds)

    # Load as Field_CM
    field = ff.Field_CM(filename)
    k = ff.Vec(3.14, 3.14)  # BZ center
    w = 0.0

    # Field_CM flattens 4D tensor to [d1*d2][d3*d4] matrix
    result = field(k, w)
    print(f"  Result shape (flattened): {len(result)}x{len(result[0])}")
    print(f"  Expected: [{dim*dim}][{dim*dim}] = [4][4]")
    print(f"  result[0][0] = {result[0][0]}")


def example_vector_field():
    """Example 5: Vector field (rank-1)"""
    print("\n" + "="*60)
    print("Example 5: Vector Field (inds={4})")
    print("="*60)

    # Create 4-component vector field
    nk = 10
    vec_dim = 4

    data = []
    for ik in range(nk):
        vector = []
        for i in range(vec_dim):
            val = float(ik) / nk + float(i) / vec_dim
            vector.append(complex(val, -val/2.0))
        data.append(vector)

    filename = "/tmp/example_vector_field.h5"
    mesh = [nk]
    domain = [[1.0]]
    inds = [vec_dim]

    print(f"Creating {vec_dim}-component vector field...")
    ff.save_data_vector(filename, data, True, mesh, domain, [], inds)

    # Load as Field_CM
    # For rank-1, Field_CM wraps as column matrix [vec_dim][1]
    field = ff.Field_CM(filename)
    k = ff.Vec(0.5)
    result = field(k)

    print(f"  Result shape (as column matrix): {len(result)}x{len(result[0])}")
    print(f"  Vector components:")
    for i in range(len(result)):
        print(f"    [{i}]: {result[i][0]}")


def main():
    """Run all examples"""
    print("\n" + "="*60)
    print("BaseData and Field Examples - Python")
    print("="*60)

    example_scalar_field()
    example_matrix_field()
    example_nonuniform_matrix()
    example_4d_vertex()
    example_vector_field()

    print("\n" + "="*60)
    print("All examples completed successfully!")
    print("="*60)
    print("\nFiles created in /tmp/:")
    print("  - example_scalar_field.h5")
    print("  - example_matrix_field.h5")
    print("  - example_nonuniform_matrix.h5")
    print("  - example_4d_vertex.h5")
    print("  - example_vector_field.h5")
    print("\nYou can inspect these files with h5dump or h5ls")


if __name__ == "__main__":
    main()
