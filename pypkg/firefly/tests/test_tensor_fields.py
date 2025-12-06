"""
Test script for tensor field save/load functionality in Python.
Tests both 3D and 4D tensor fields.
"""
import numpy as np
import sys
import os

# Add firefly to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.module.imports.cpp_imports import (
    BaseData, save_data_tensor3, save_data_tensor4
)

def test_tensor3_save_load():
    """Test 3D tensor (n_indices=3) save and load."""
    print("Testing 3D tensor save/load...")

    # Create test data: 3 spatial points, 2x2x2 tensor at each point
    num_tensors = 3
    ten_dim = 2
    mesh = np.array([3], dtype=np.int32)
    domain = np.array([[1.0]], dtype=np.float32)
    w_points = np.array([], dtype=np.float32)

    # Create tensor data: shape (num_tensors, dim, dim, dim)
    data = np.zeros((num_tensors, ten_dim, ten_dim, ten_dim), dtype=np.complex64)
    for t in range(num_tensors):
        for i in range(ten_dim):
            for j in range(ten_dim):
                for k in range(ten_dim):
                    val = float(t + i + j + k)
                    data[t, i, j, k] = complex(val, val / 10.0)

    # Save
    filename = "/tmp/test_tensor3_python.h5"
    save_data_tensor3(filename, data, num_tensors, ten_dim, True, mesh, domain, w_points)
    print(f"  Saved 3D tensor to {filename}")

    # Load
    loaded = BaseData(filename)
    print(f"  Loaded metadata: n_indices={loaded.n_indices}, dim_indices={loaded.dim_indices}")

    # Check metadata
    assert loaded.n_indices == 3, f"Expected n_indices=3, got {loaded.n_indices}"
    assert loaded.dim_indices == ten_dim, f"Expected dim_indices={ten_dim}, got {loaded.dim_indices}"
    assert loaded.nk == num_tensors, f"Expected nk={num_tensors}, got {loaded.nk}"

    # Get data back
    loaded_data = loaded.get_data()
    print(f"  Loaded data shape: {loaded_data.shape}")

    # Check shape
    expected_shape = (num_tensors, ten_dim, ten_dim, ten_dim)
    assert loaded_data.shape == expected_shape, f"Expected shape {expected_shape}, got {loaded_data.shape}"

    # Check values
    max_diff = np.max(np.abs(loaded_data - data))
    print(f"  Max difference: {max_diff}")
    assert max_diff < 1e-5, f"Data mismatch, max diff = {max_diff}"

    # Cleanup
    os.remove(filename)
    print("  ✓ 3D tensor test passed!")
    return True

def test_tensor4_save_load():
    """Test 4D tensor (n_indices=4) save and load."""
    print("Testing 4D tensor save/load...")

    # Create test data: 4 spatial points, 2x2x2x2 tensor at each point
    num_tensors = 4
    ten_dim = 2
    mesh = np.array([2, 2], dtype=np.int32)
    domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    w_points = np.array([], dtype=np.float32)

    # Create tensor data: shape (num_tensors, dim, dim, dim, dim)
    data = np.zeros((num_tensors, ten_dim, ten_dim, ten_dim, ten_dim), dtype=np.complex64)
    for t in range(num_tensors):
        for i in range(ten_dim):
            for j in range(ten_dim):
                for k in range(ten_dim):
                    for l in range(ten_dim):
                        val = float(t + i + j + k + l)
                        data[t, i, j, k, l] = complex(val, val / 10.0)

    # Save
    filename = "/tmp/test_tensor4_python.h5"
    save_data_tensor4(filename, data, num_tensors, ten_dim, True, mesh, domain, w_points)
    print(f"  Saved 4D tensor to {filename}")

    # Load
    loaded = BaseData(filename)
    print(f"  Loaded metadata: n_indices={loaded.n_indices}, dim_indices={loaded.dim_indices}")

    # Check metadata
    assert loaded.n_indices == 4, f"Expected n_indices=4, got {loaded.n_indices}"
    assert loaded.dim_indices == ten_dim, f"Expected dim_indices={ten_dim}, got {loaded.dim_indices}"
    assert loaded.nk == num_tensors, f"Expected nk={num_tensors}, got {loaded.nk}"

    # Get data back
    loaded_data = loaded.get_data()
    print(f"  Loaded data shape: {loaded_data.shape}")

    # Check shape
    expected_shape = (num_tensors, ten_dim, ten_dim, ten_dim, ten_dim)
    assert loaded_data.shape == expected_shape, f"Expected shape {expected_shape}, got {loaded_data.shape}"

    # Check values
    max_diff = np.max(np.abs(loaded_data - data))
    print(f"  Max difference: {max_diff}")
    assert max_diff < 1e-5, f"Data mismatch, max diff = {max_diff}"

    # Cleanup
    os.remove(filename)
    print("  ✓ 4D tensor test passed!")
    return True

def test_tensor4_with_frequency():
    """Test 4D tensor with frequency dimension."""
    print("Testing 4D tensor with frequency dimension...")

    # Create test data: 2 w-points × 3 k-points, 2x2x2x2 tensor at each
    nw = 2
    nk = 3
    num_tensors = nw * nk
    ten_dim = 2
    mesh = np.array([3], dtype=np.int32)
    domain = np.array([[1.0]], dtype=np.float32)
    w_points = np.array([0.0, 1.0], dtype=np.float32)

    # Create tensor data: shape (num_tensors, dim, dim, dim, dim)
    # In k-w ordering (C++ default): [k0w0, k0w1, k1w0, k1w1, k2w0, k2w1]
    data = np.zeros((num_tensors, ten_dim, ten_dim, ten_dim, ten_dim), dtype=np.complex64)
    for k in range(nk):
        for w in range(nw):
            t = k * nw + w
            for i in range(ten_dim):
                for j in range(ten_dim):
                    for kk in range(ten_dim):
                        for l in range(ten_dim):
                            val = float(w + k + i + j + kk + l)
                            data[t, i, j, kk, l] = complex(val, val / 10.0)

    # Save (defaults to k-w ordering)
    filename = "/tmp/test_tensor4_w_python.h5"
    save_data_tensor4(filename, data, num_tensors, ten_dim, True, mesh, domain, w_points)
    print(f"  Saved 4D tensor with frequency to {filename}")

    # Load with k-w ordering (default)
    loaded = BaseData(filename, ordering="k-w")
    print(f"  Loaded metadata: n_indices={loaded.n_indices}, nw={loaded.nw}, nk={loaded.nk}")

    # Check metadata
    assert loaded.n_indices == 4, f"Expected n_indices=4, got {loaded.n_indices}"
    assert loaded.nw == nw, f"Expected nw={nw}, got {loaded.nw}"
    assert loaded.nk == nk, f"Expected nk={nk}, got {loaded.nk}"

    # Get data back
    loaded_data = loaded.get_data()
    print(f"  Loaded data shape: {loaded_data.shape}")

    # Check values
    max_diff = np.max(np.abs(loaded_data - data))
    print(f"  Max difference: {max_diff}")
    assert max_diff < 1e-5, f"Data mismatch, max diff = {max_diff}"

    # Cleanup
    os.remove(filename)
    print("  ✓ 4D tensor with frequency test passed!")
    return True

if __name__ == "__main__":
    print("=" * 60)
    print("Tensor Field Python Interface Tests")
    print("=" * 60)

    try:
        test_tensor3_save_load()
        print()
        test_tensor4_save_load()
        print()
        test_tensor4_with_frequency()
        print()
        print("=" * 60)
        print("✓ All Python tensor field tests passed!")
        print("=" * 60)
        sys.exit(0)
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
