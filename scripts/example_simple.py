#!/usr/bin/env python3
"""
BaseData and Fields Tutorial. Very simple, just go to bottom of file
"""

import os
import numpy as np
import firefly as fly

file = "temp_data.h5"
# f(w,k) = exp(-k^2/(0.05+w)^2)
def gaussian(k, w):
    kmag = k[0]**2 + k[1]**2
    sigma = 0.2 + w
    return np.exp(-kmag / sigma**2)

def get_data():
    # w-points are defined, for this example k-points are assumed to lie on a mesh
    # If non-mesh functions are desired, they can also be passed in as a list of points
    w_points = [0.0, 0.5, 1.0]
    nw = len(w_points)
    nk = 100
    k_mesh = [nk, nk] # 2D Mesh
    print("w-points:", w_points)
    print("k-mesh:", k_mesh)

    # K-space mesh
    kx = np.linspace(-0.5, 0.5, nk)
    ky = np.linspace(-0.5, 0.5, nk)
    kgrid = np.meshgrid(kx, ky, indexing='ij')
    kgrid = np.stack(kgrid, axis=-1).reshape(-1, 2)  # Shape (nk*nk, 2)
    # Brillouin Zone (in 3D by default)
    BZ = np.array([[2*np.pi, 0.0, 0.0], [0.0, 2*np.pi, 0.0], [0.0, 0.0, 2*np.pi]]) 
    # Brillouin Zone changed to 2D
    BZ = BZ[:2, :2]
    kpts = kgrid @ BZ
    print("kmin:", kpts.min(axis=0))
    print("kmax:", kpts.max(axis=0))
    print("K-point list shape:", kpts.shape)

    # Create Gaussian data: f(w,k) = exp(-k^2/(0.05+w)^2)
    # w-k ordering: frequency varies slowest
    print(f"\nCreating Gaussian data on {nw}x{nk} grid...")
    #data = np.zeros((nw, nk, nk), dtype=complex)
    data = np.zeros((nw, nk, nk))
    for i in range(nw):
        data[i, :, :] = gaussian(kpts.T, w_points[i]).reshape(nk, nk)
    print("data shape:", data.shape)
    return data, k_mesh, BZ, w_points

def BaseData_example(data, k_mesh, BZ, w_points):
    # Step 1: Save using save_data_scalar
    print("\nStep 1: Save data")
    fly.save_data(file, data, mesh=k_mesh, domain=BZ, w_points=w_points)
    #fly.save_data_scalar(file1, data, True, mesh, domain, w_points)
    print(f"  Saved to: {file}")

    # Step 2: Load into BaseData
    print("\nStep 2: Load into BaseData")
    bd = fly.BaseData(file)
    print(f"  Loaded BaseData:")
    print(f"    is_complex: {bd.is_complex}")
    print(f"    inds: {bd.inds}")
    print(f"    rank: {len(bd.inds)}")
    print(f"    mesh: {bd.mesh}")
    print(f"    nk: {bd.nk}")
    print(f"    nw: {bd.nw}")

    # Access data through the data property
    print(f"\n  Accessing data through bd.data property:")
    print(f"    data shape: {bd.data.shape}")
    print(f"    data dtype: {bd.data.dtype}")
    print(f"    data min/max: {bd.data.min():.6f} / {bd.data.max():.6f}")

    # Step 3: Save BaseData
    print("\nStep 3: Save BaseData")
    bd.save(file)
    return bd

def Field_R_example():
    print("\nLoad into Field_R (Real Scalar Field) from file")
    field = fly.Field_R(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]  
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print(f"    Evaluating at k=(0.0,0.0), w=0.5:", end='')
    print(f"    Expected, found: {expected:.6f}, {val}")

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5] 
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print(f"    Evaluating at k=(-0.2, 0.5), w=1.0:", end='')
    print(f"    Expected, found: {expected}, {val}")

def Field_C_example():
    # Step 4: Load into Field_R (Real Scalar Field) from file
    print("\nLoad into Field_C (Complex Scalar Field) from file")
    field = fly.Field_C(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]  
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print(f"    Evaluating at k=(0.0,0.0), w=0.5:", end='')
    print(f"    Expected, found: {expected}, {val}")

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5] 
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print(f"    Evaluating at k=(-0.2, 0.5), w=1.0:", end='')
    print(f"    Expected, found: {expected}, {val}")

def Field_CM_example():
    # Step 4: Load into Field_R (Real Scalar Field) from file
    print("\nLoad into Field_CM (Complex Matrix Field) from file")
    field = fly.Field_CM(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]  
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print(f"    Evaluating at k=(0.0,0.0), w=0.5:", end='')
    print(f"    Expected, found: {expected}, {val}")

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5] 
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print(f"    Evaluating at k=(-0.2, 0.5), w=1.0:", end='')
    print(f"    Expected, found: {expected}, {val}")

def Field_RM_example():
    # Step 4: Load into Field_R (Real Scalar Field) from file
    print("\nLoad into Field_RM (Complex Matrix Field) from file")
    field = fly.Field_RM(file)

    # Test at k=(0,0), w=0.5
    k_test = [0.0, 0.0]  
    w_test = 0.5

    val = field(k_test, w_test)
    expected = gaussian(k_test, w_test)

    print(f"    Evaluating at k=(0.0,0.0), w=0.5:", end='')
    print(f"    Expected, found: {expected}, {val}")

    # Test at k=(-0.2,0.5), w=1.0
    k_test = [-0.2, 0.5] 
    w_test = 1.0
    expected = gaussian(k_test, w_test)

    val = field(k_test, w_test)

    print(f"    Evaluating at k=(-0.2, 0.5), w=1.0:", end='')
    print(f"    Expected, found: {expected}, {val}")

if __name__ == "__main__":
    # Create data on a 3D w-k grid
    data, k_mesh, BZ, w_points = get_data() 
    fly.save_data(file, data, mesh=k_mesh, domain=BZ, w_points=w_points) # Save data

    bd = BaseData_example(data, k_mesh, BZ, w_points) # Load into BaseData
    bd.save(file) # Save again

    # Make a Real Scalar Field
    Field_R_example()

    # Make Data Complex
    bd.data = bd.data.astype(complex)
    bd.save(file)
    # Make a Complex Scalar Field
    Field_C_example()

    bd.data = bd.data.reshape(bd.nw, bd.mesh[0], bd.mesh[1], 1, 1)
    bd.inds = [1,1]
    bd.save(file)

    Field_CM_example()

    bd.data = bd.data.real.astype(float)
    bd.data = bd.data.reshape(bd.nw, bd.mesh[0], bd.mesh[1], 1, 1)
    bd.save(file)

    Field_RM_example()

    if os.path.exists(file):
        os.remove(file)
        print(f"\nFile '{file}' deleted successfully.")

