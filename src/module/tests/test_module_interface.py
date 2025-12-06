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

# Round-trip tests - save and read back
def test_save_read_scalar_real():
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_roundtrip_scalar_real.h5")

        # Create test data - 10x10 grid with values = x + y
        nk1, nk2 = 10, 10
        data = np.zeros((nk1, nk2), dtype=np.float32)
        for i in range(nk1):
            for j in range(nk2):
                data[i, j] = float(i + j)

        mesh = [nk1, nk2]
        domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)

        # Save the data
        ff.save_data_scalar(filename, data, False, mesh, domain)

        # Read it back using Field_R
        field = ff.Field_R(filename)

        # Verify data at interior points only (avoid boundaries)
        passed = True
        tolerance = 1e-3
        for i in [1, 3, 5, 7]:
            for j in [1, 3, 5, 7]:
                k_point = [float(i)/(nk1-1) - 0.5, float(j)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(value - expected) > tolerance:
                    print(f"Mismatch at ({i},{j}): got {value}, expected {expected}")
                    passed = False

        import shutil
        shutil.rmtree(tmpdir)
        return passed
    except Exception as e:
        print(f"save_read_scalar_real error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_save_read_scalar_complex():
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_roundtrip_scalar_complex.h5")

        # Create test data - 8x8 grid with complex values
        nk1, nk2 = 8, 8
        data = np.zeros((nk1, nk2), dtype=np.complex64)
        for i in range(nk1):
            for j in range(nk2):
                data[i, j] = complex(float(i), float(j))

        mesh = [nk1, nk2]
        domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)

        # Save the data
        ff.save_data_scalar(filename, data, True, mesh, domain)

        # Read it back using Field_C
        field = ff.Field_C(filename)

        # Verify data at interior points only
        passed = True
        tolerance = 1e-3
        for i in [1, 3, 5]:
            for j in [1, 3, 5]:
                k_point = [float(i)/(nk1-1) - 0.5, float(j)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(value.real - expected.real) > tolerance or abs(value.imag - expected.imag) > tolerance:
                    print(f"Mismatch at ({i},{j}): got {value}, expected {expected}")
                    passed = False

        import shutil
        shutil.rmtree(tmpdir)
        return passed
    except Exception as e:
        print(f"save_read_scalar_complex error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_save_read_with_frequency():
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_roundtrip_freq.h5")

        # Create test data - 5x5 k-grid with 8 frequencies
        # FFirefly convention: (nw, nx, ny)
        nk1, nk2 = 5, 5
        nw = 8
        data = np.zeros((nw, nk1, nk2), dtype=np.float32)
        for w in range(nw):
            for i in range(nk1):
                for j in range(nk2):
                    data[w, i, j] = float(i + j + w)

        mesh = [nk1, nk2]  # Only k-space dimensions
        domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
        w_points = np.array([float(w) for w in range(nw)], dtype=np.float32)

        # Save the data
        ff.save_data_scalar(filename, data, False, mesh, domain, w_points)

        # Read it back using Field_R
        field = ff.Field_R(filename)

        # Verify data at interior points with different frequencies
        passed = True
        tolerance = 1e-3
        for i in [1, 2, 3]:
            for j in [1, 2, 3]:
                for w in [1, 3, 5]:
                    k_point = [float(i)/(nk1-1) - 0.5, float(j)/(nk2-1) - 0.5, 0.0]
                    value = field(k_point, float(w))
                    expected = data[w, i, j]  # Changed from data[i,j,w] to match (nw,nx,ny) layout
                    if abs(value - expected) > tolerance:
                        print(f"Mismatch at ({i},{j},w={w}): got {value}, expected {expected}")
                        passed = False

        import shutil
        shutil.rmtree(tmpdir)
        return passed
    except Exception as e:
        print(f"save_read_with_frequency error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_save_data_dispatcher_real():
    """Test save_data() automatically dispatches for real data"""
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_dispatcher_real.h5")

        # Create simple real data
        nk1, nk2 = 6, 6
        data = np.zeros((nk1, nk2), dtype=np.float32)
        for i in range(nk1):
            for j in range(nk2):
                data[i, j] = float(i + j)

        mesh = [nk1, nk2]
        domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)

        # Use save_data (dispatcher) instead of save_data_scalar
        ff.save_data(filename, data, mesh=mesh, domain=domain)

        # Read back and verify (interior points only)
        field = ff.Field_R(filename)
        passed = True
        tolerance = 1e-3
        for i in [1, 2, 3, 4]:
            for j in [1, 2, 3, 4]:
                k_point = [float(i)/(nk1-1) - 0.5, float(j)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(value - expected) > tolerance:
                    print(f"Dispatcher real mismatch at ({i},{j}): got {value}, expected {expected}")
                    passed = False

        import shutil
        shutil.rmtree(tmpdir)
        return passed
    except Exception as e:
        print(f"save_data_dispatcher_real error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_save_data_dispatcher_complex():
    """Test save_data() automatically dispatches for complex data"""
    try:
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_dispatcher_complex.h5")

        # Create complex data
        nk1, nk2 = 6, 6
        data = np.zeros((nk1, nk2), dtype=np.complex64)
        for i in range(nk1):
            for j in range(nk2):
                data[i, j] = complex(float(i), float(j))

        mesh = [nk1, nk2]
        domain = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)

        # Use save_data (dispatcher) - should auto-detect complex
        ff.save_data(filename, data, mesh=mesh, domain=domain)

        # Read back and verify (interior points only)
        field = ff.Field_C(filename)
        passed = True
        tolerance = 1e-3
        for i in [1, 2, 3, 4]:
            for j in [1, 2, 3, 4]:
                k_point = [float(i)/(nk1-1) - 0.5, float(j)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(value.real - expected.real) > tolerance or abs(value.imag - expected.imag) > tolerance:
                    print(f"Dispatcher complex mismatch at ({i},{j}): got {value}, expected {expected}")
                    passed = False

        import shutil
        shutil.rmtree(tmpdir)
        return passed
    except Exception as e:
        print(f"save_data_dispatcher_complex error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_unified_field_scalar_real():
    """Test Field with real scalar data"""
    try:
        import tempfile
        import os
        import shutil
        
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_unified_scalar_real.h5")
        
        # Create real scalar data using save_data
        nk = 6
        data = np.array([float(i) for i in range(nk)], dtype=np.float32)
        mesh = [nk]
        domain = np.array([[1.0]], dtype=np.float32)
        ff.save_data(filename, data, mesh=mesh, domain=domain)

        # Load with Field
        field = ff.Field(filename)

        # Check type flags
        if field.is_complex or field.is_matrix:
            print(f"Field type flags incorrect: is_complex={field.is_complex}, is_matrix={field.is_matrix}")
            shutil.rmtree(tmpdir)
            return False

        # Test evaluation
        value = field([0.5, 0.0, 0.0])
        if not isinstance(value, (float, np.floating)):
            print(f"Field scalar real returned wrong type: {type(value)}")
            shutil.rmtree(tmpdir)
            return False
        
        shutil.rmtree(tmpdir)
        return True
    except Exception as e:
        print(f"test_unified_field_scalar_real error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_unified_field_scalar_complex():
    """Test Field with complex scalar data"""
    try:
        import tempfile
        import os
        import shutil
        
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_unified_scalar_complex.h5")
        
        # Create complex scalar data
        nk = 6
        data = np.array([complex(i, i+1) for i in range(nk)], dtype=np.complex64)
        mesh = [nk]
        domain = np.array([[1.0]], dtype=np.float32)
        ff.save_data(filename, data, mesh=mesh, domain=domain)

        # Load with Field
        field = ff.Field(filename)

        # Check type flags
        if not field.is_complex or field.is_matrix:
            print(f"Field type flags incorrect: is_complex={field.is_complex}, is_matrix={field.is_matrix}")
            shutil.rmtree(tmpdir)
            return False

        # Test evaluation
        value = field([0.5, 0.0, 0.0])
        if not isinstance(value, (complex, np.complexfloating)):
            print(f"Field scalar complex returned wrong type: {type(value)}")
            shutil.rmtree(tmpdir)
            return False
        
        shutil.rmtree(tmpdir)
        return True
    except Exception as e:
        print(f"test_unified_field_scalar_complex error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_unified_field_matrix_real():
    """Test Field with real matrix data"""
    try:
        import tempfile
        import os
        import shutil
        
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_unified_matrix_real.h5")
        
        # Create real matrix data (2x2 matrices at each point)
        nk = 4
        nbnd = 2
        data = np.zeros((nk, nbnd, nbnd), dtype=np.float32)
        for i in range(nk):
            data[i] = np.eye(nbnd, dtype=np.float32) * (i + 1)

        mesh = [nk]
        domain = np.array([[1.0]], dtype=np.float32)
        ff.save_data(filename, data, mesh=mesh, domain=domain, dim_indices=nbnd)

        # Load with Field
        field = ff.Field(filename)

        # Check type flags
        if field.is_complex or not field.is_matrix:
            print(f"Field type flags incorrect: is_complex={field.is_complex}, is_matrix={field.is_matrix}")
            shutil.rmtree(tmpdir)
            return False

        # Test evaluation
        value = field([0.5, 0.0, 0.0])
        if not isinstance(value, np.ndarray) or value.dtype != np.float32:
            print(f"Field matrix real returned wrong type: {type(value)}, dtype={value.dtype if hasattr(value, 'dtype') else 'N/A'}")
            shutil.rmtree(tmpdir)
            return False

        if value.shape != (nbnd, nbnd):
            print(f"Field matrix real returned wrong shape: {value.shape}")
            shutil.rmtree(tmpdir)
            return False
        
        shutil.rmtree(tmpdir)
        return True
    except Exception as e:
        print(f"test_unified_field_matrix_real error: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_unified_field_matrix_complex():
    """Test Field with complex matrix data"""
    try:
        import tempfile
        import os
        import shutil
        
        tmpdir = tempfile.mkdtemp()
        filename = os.path.join(tmpdir, "test_unified_matrix_complex.h5")
        
        # Create complex matrix data (2x2 matrices at each point)
        nk = 4
        nbnd = 2
        data = np.zeros((nk, nbnd, nbnd), dtype=np.complex64)
        for i in range(nk):
            data[i] = np.eye(nbnd, dtype=np.complex64) * complex(i + 1, i + 2)

        mesh = [nk]
        domain = np.array([[1.0]], dtype=np.float32)
        ff.save_data(filename, data, mesh=mesh, domain=domain, dim_indices=nbnd)

        # Load with Field
        field = ff.Field(filename)

        # Check type flags
        if not field.is_complex or not field.is_matrix:
            print(f"Field type flags incorrect: is_complex={field.is_complex}, is_matrix={field.is_matrix}")
            shutil.rmtree(tmpdir)
            return False

        # Test evaluation
        value = field([0.5, 0.0, 0.0])
        if not isinstance(value, np.ndarray) or value.dtype != np.complex64:
            print(f"Field matrix complex returned wrong type: {type(value)}, dtype={value.dtype if hasattr(value, 'dtype') else 'N/A'}")
            shutil.rmtree(tmpdir)
            return False

        if value.shape != (nbnd, nbnd):
            print(f"Field matrix complex returned wrong shape: {value.shape}")
            shutil.rmtree(tmpdir)
            return False
        
        shutil.rmtree(tmpdir)
        return True
    except Exception as e:
        print(f"test_unified_field_matrix_complex error: {e}")
        import traceback
        traceback.print_exc()
        return False
