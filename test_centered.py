#!/usr/bin/env python3
"""Test script for centered parameter in save_data functions."""

import sys
import os
import numpy as np

# Add the module path
sys.path.insert(0, '/home/g/Research/FFirefly/src/module/imports')
from cpp_imports import save_data, BaseData

def test_centered():
    """Test that centered parameter is saved and loaded correctly."""
    # Create test data
    mesh = np.array([4, 4], dtype=np.int32)
    domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    data = np.random.rand(4, 4).astype(np.float32)

    # Test 1: Save with centered=True (default)
    save_data("/tmp/test_centered_true.h5", data, mesh=mesh, domain=domain, centered=True)
    bd_true = BaseData("/tmp/test_centered_true.h5")
    print(f"Test 1 - centered=True: loaded centered = {bd_true.centered}")
    assert bd_true.centered == True, "Expected centered=True"

    # Test 2: Save with centered=False
    save_data("/tmp/test_centered_false.h5", data, mesh=mesh, domain=domain, centered=False)
    bd_false = BaseData("/tmp/test_centered_false.h5")
    print(f"Test 2 - centered=False: loaded centered = {bd_false.centered}")
    assert bd_false.centered == False, "Expected centered=False"

    # Test 3: Default value (should be True)
    save_data("/tmp/test_centered_default.h5", data, mesh=mesh, domain=domain)
    bd_default = BaseData("/tmp/test_centered_default.h5")
    print(f"Test 3 - default: loaded centered = {bd_default.centered}")
    assert bd_default.centered == True, "Expected default centered=True"

    # Cleanup
    for f in ["/tmp/test_centered_true.h5", "/tmp/test_centered_false.h5", "/tmp/test_centered_default.h5"]:
        if os.path.exists(f):
            os.remove(f)

    print("\nAll Python tests passed!")

if __name__ == "__main__":
    test_centered()
