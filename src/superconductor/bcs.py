import firefly as fly
import firefly.config as cfg
from firefly.diagram import *

from triqs.gf import *
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.lattice import BrillouinZone, BravaisLattice
from triqs.gf.mesh_product import MeshProduct
import numpy as np
from triqs_tprf.lattice import *
from triqs_tprf import *

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

def bcs():
    H_r, kmesh, e_k = fly.load_triqs_H.get_energy_mesh()
    Delta = Gf(name='Delta', mesh=kmesh, target_shape=[nstates, nstates])
    V = Gf(name='Vertex', mesh=kmesh, target_shape=[nstates, nstates, nstates, nstates])
    kmesh = get_k_mesh(BZ, nx, ny, nz)
    print(outdir + prefix + '_vertex.h5')
    vertex = fly.Field_CM(outdir + prefix + '_vertex.h5')
    print(vertex([0.1,0.2,0.3]))
    data = vertex(kmesh)
    V.data[:] = np.reshape(data, (Nk, nstates, nstates, nstates, nstates))

    # Randomize initial gap function
    #Delta.data[:] = 1/(nx*ny*nz) * (np.random.rand(Nk, nstates, nstates) - 0.5)
