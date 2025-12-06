import firefly as fly
import firefly.config as cfg
from firefly.diagram import *

from triqs.gf import *
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.lattice import BrillouinZone, BravaisLattice
from triqs.gf.mesh_product import MeshProduct
from triqs_tprf.lattice import *
from triqs_tprf import *

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


outdir = cfg.outdir
prefix = cfg.prefix

nstates = cfg.nstates
Nk = cfg.k_mesh[0]
nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim == 2:
    nz = 1
BZ = np.array(cfg.brillouin_zone)

mu = cfg.fermi_energy
beta = 1.0 / cfg.Temperature

def get_k_mesh(BZ, nx, ny, nz):
    # fractional grid (0→1)
    fx = np.linspace(0, 1, nx, endpoint=False)
    fy = np.linspace(0, 1, ny, endpoint=False)
    fz = np.linspace(0, 1, nz, endpoint=False)

    # 3D fractional mesh
    F = np.stack(np.meshgrid(fx, fy, fz, indexing='ij'), axis=-1)  # (nx,ny,nz,3)

    # convert to Cartesian k = f1*b1 + f2*b2 + f3*b3
    K = F @ BZ.T                          # shape (nx,ny,nz,3)

    return K.reshape(-1, 3)           # full grid, flattened grid

def make_lanczos(A):

    def mv(x):
        #y = np.zeros_like(x)
        y = np.fft(x)
        convd = A * y
        result = np.ifft(convd)
        return result

    n = 1000
    A = LinearOperator((n, n), matvec=mv, dtype=float)

    # ARPACK uses Lanczos for symmetric problems
    vals, vecs = eigsh(A, k=5, which="LM")  # 5 largest-magnitude eigenpairs

def bcs():
    H_r, kmesh, e_k = fly.load_triqs_H.get_energy_mesh()
    Delta = Gf(name='Delta', mesh=kmesh, target_shape=[nstates, nstates])
    V = Gf(name='Vertex', mesh=kmesh, target_shape=[nstates, nstates, nstates, nstates])
    kmesh = get_k_mesh(BZ, nx, ny, nz)

    # Load vertex data using BaseData (4-index tensor stored as scalar for single-band)
    vertex_file = outdir + prefix + '_vertex.h5'
    print(f"Loading vertex from {vertex_file}")
    test_V = fly.Field_CM(vertex_file)
    print(test_V([0.0,0.0,0.0]))
    base_data = fly.BaseData(vertex_file)

    # Get metadata
    n_indices = base_data.n_indices
    dim_indices = base_data.dim_indices
    print(f"Vertex file: n_indices={n_indices}, dim_indices={dim_indices}")

    # Get raw data
    data_flat = base_data.get_data()
    nw = len(base_data.w_points)
    nk_total = len(data_flat) // nw
    mesh_shape = tuple(base_data.mesh)
    nk = np.prod(mesh_shape)  # Total number of k-points
    print(f"Data shape: nw={nw}, nk={nk}, mesh={mesh_shape}")

    # Reshape to (nw, nk) for single-band system
    vertex_data = np.array(data_flat).reshape((nw, nk))

    # Take w=0 component (middle frequency for now, should find closest to 0)
    w0_idx = nw // 2
    V.data[:, 0, 0, 0, 0] = vertex_data[w0_idx, :]

    eigs, vecs = make_lanczos(V.data[:,0,0,0,0])
    print(eigs)
    for i in range(len(eigs)):
        if eigs[i] > 0:
            fly.save_data(outdir + prefix + "_gap.h5", vecs[i], cfg.k_mesh, BZ)
            break
    print("Saved to ", outdir + prefix + "_gap.h5")
