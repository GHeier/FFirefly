"""Test module for C++ exports to Python interface - Basic functionality tests."""

import sys
sys.path.insert(0, "/home/g/Research/FFirefly/src/module/imports")

import cpp_imports as ff
import numpy as np
import os
import tempfile

# Vec class tests
def test_vec_constructor_empty():
    try:
        v = ff.Vec()
        return True
    except Exception as e:
        print(f"Vec empty constructor error: {e}")
        return False

def test_vec_constructor_args():
    try:
        v = ff.Vec(1.0, 2.0, 3.0)
        result = (abs(v.x - 1.0) < 0.001 and abs(v.y - 2.0) < 0.001 and abs(v.z - 3.0) < 0.001)
        return result
    except Exception as e:
        print(f"Vec args constructor error: {e}")
        return False

def test_vec_getters():
    try:
        v = ff.Vec(1.5, 2.5, 3.5, 0.0, 0.0, 3, 0)
        return (abs(v.x - 1.5) < 0.001 and
                abs(v.y - 2.5) < 0.001 and
                abs(v.z - 3.5) < 0.001 and
                v.dimension == 3)
    except Exception as e:
        print(f"Vec getters error: {e}")
        return False

# Load config test
def test_load_config():
    try:
        config_path = "/home/g/Research/FFirefly/build/bin/input.cfg"
        if os.path.exists(config_path):
            ff.load_config(config_path)
            return True
        return False
    except Exception as e:
        print(f"load_config error: {e}")
        return False

# Save data tests
def test_save_data_scalar():
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_scalar.h5")
        data = np.ones((10, 10), dtype=np.float32)
        mesh = [10, 10]
        domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
        ff.save_data_scalar(filename, data, False, mesh, domain)
        result = os.path.exists(filename)
        import shutil
        shutil.rmtree(tmpdir)
        return result
    except Exception as e:
        print(f"save_data_scalar error: {e}")
        return False

def test_save_data_scalar_complex():
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_scalar_complex.h5")
        data = np.ones((10, 10), dtype=np.complex64)
        mesh = [10, 10]
        domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
        ff.save_data_scalar(filename, data, True, mesh, domain)
        result = os.path.exists(filename)
        import shutil
        shutil.rmtree(tmpdir)
        return result
    except Exception as e:
        print(f"save_data_scalar complex error: {e}")
        return False

def test_save_data_vector():
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_vector.h5")
        nk = 10
        vec_len = 3
        data = np.ones((nk, vec_len), dtype=np.float32)
        mesh = [10]
        domain = np.array([[1.0]], dtype=np.float32)
        ff.save_data_vector(filename, data, nk, vec_len, False, mesh, domain)
        result = os.path.exists(filename)
        import shutil
        shutil.rmtree(tmpdir)
        return result
    except Exception as e:
        print(f"save_data_vector error: {e}")
        return False

def test_save_data_matrix():
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_matrix.h5")
        mat_dim = 2
        num_matrices = 5
        data = np.ones((num_matrices, mat_dim, mat_dim), dtype=np.float32)
        mesh = [5]
        domain = np.array([[1.0]], dtype=np.float32)
        ff.save_data_matrix(filename, data, num_matrices, mat_dim, False, mesh, domain)
        result = os.path.exists(filename)
        import shutil
        shutil.rmtree(tmpdir)
        return result
    except Exception as e:
        print(f"save_data_matrix error: {e}")
        return False
