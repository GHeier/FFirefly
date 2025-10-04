import numpy as np
import tbmodels as tb
import itertools
import matplotlib.pyplot as plt

import firefly as fly
import firefly.config as cfg

outdir = cfg.outdir
prefix = cfg.prefix

nstates = cfg.nstates
nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
BZ = cfg.brillouin_zone
if dim == 2:
    nz = 1
    BZ = BZ[:2, :2]

t0 = cfg.t0



# Create Mesh
kx = (np.arange(nx) + 0.5) / nx
ky = (np.arange(ny) + 0.5) / ny
kz = (np.arange(nz) + 0.5) / nz
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
kpts = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], axis=1)




# TIGHT BINDING 

if cfg.hamiltonian == "tight_binding":
    # one orbital per cell at the origin
    model = tb.Model(on_site=[0.0], dim=dim, pos=[[0.0]*dim])

    # H(R) entries for NN on a square lattice, basis size = 1
    # add both R and -R to keep H Hermitian
    R = [(1,0,0), (0,1,0), (0,0,1)]
    for i in range(dim):
        model.add_hop(-t0, 0, 0, R[i])

    Hk = model.hamilton(k=kpts).reshape((nx, ny, nz, nstates, nstates))
    fly.save_data(outdir + prefix + "_Hk.h5", Hk, True, [nx, ny, nz], BZ)


#
## optional: band on a simple path Γ→X→M→Γ
#def kline(a, b, n):
#    a, b = np.array(a), np.array(b)
#    t = np.linspace(0, 1, n, endpoint=False)
#    return [list((1 - s) * a + s * b) for s in t]
#
## symmetry points
#G = [0.0, 0.0, 0.0]
#X = [0.5, 0.0, 0.0]
#M = [0.5, 0.5, 0.0]
#
## path segments
#path = kline(G, X, 100) + kline(X, M, 100) + kline(M, G, 100) + [G]
#
## compute band energies
#bands = np.array([model.eigenval(k=kpt) for kpt in path]).squeeze()
#
## build 1D coordinate along the path
#k_coords = np.cumsum([0] + [np.linalg.norm(np.array(path[i+1]) - np.array(path[i]))
#                            for i in range(len(path)-1)])
#
#plt.figure(figsize=(5,4))
#plt.plot(k_coords, bands, 'b-')
#plt.xticks([k_coords[0], k_coords[99], k_coords[199], k_coords[-1]], ['Γ', 'X', 'M', 'Γ'])
#plt.ylabel('Energy')
#plt.xlabel('k-path')
#plt.title('2D Square-Lattice Band')
#plt.grid(True, ls='--', alpha=0.4)
#plt.tight_layout()
#plt.show()

