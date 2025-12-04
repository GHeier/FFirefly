import firefly as fly
import firefly.config as cfg

import numpy as np

outdir = cfg.outdir
prefix = cfg.prefix

nx, ny, nz = cfg.k_mesh
if cfg.dimension == 2:
    nz = 1
ne = cfg.w_pts

BZ = np.array(cfg.brillouin_zone)

def get_e_mesh():
    band = fly.Bands()
    kx, ky, kz = np.meshgrid(np.arange(nx) / nx,
                             np.arange(ny) / ny,
                             np.arange(nz) / nz,
                             indexing='ij')
    kpts = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], axis=-1)
    kpts = kpts @ BZ.T

    e_mesh = band(kpts).ravel()
    return e_mesh

def get_DOS_gaussian():
    e_mesh = get_e_mesh()
    emin = np.min(e_mesh)
    emax = np.max(e_mesh)
    grid = np.linspace(emin, emax, ne)
    X = grid[:, None] - e_mesh[None, :]
    dos = np.sum(np.exp(-X**2 / (2 * cfg.smearing**2)), axis=1)

    area = np.trapz(dos, grid)
    if area > 0:
        dos /= area
    else:
        print("DOS area is zero or negative.")
    fly.save_data(outdir + prefix + '_DOS.h5', dos, [], [[]], grid)
    print("Saved to ", outdir + prefix + "_DOS.h5")

