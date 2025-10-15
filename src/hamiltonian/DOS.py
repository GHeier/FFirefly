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

def get_DOS_gaussian(e_mesh):
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
    return grid, dos

def get_DOS():
    band = fly.Bands()
    kx, ky, kz = np.meshgrid(np.arange(nx) / nx,
                             np.arange(ny) / ny,
                             np.arange(nz) / nz,
                             indexing='ij')
    kpts = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], axis=-1)
    kpts = kpts @ BZ.T

    e_mesh = band(kpts).ravel()

    if cfg.method == 'gaussian':
        w_pts, dos = get_DOS_gaussian(e_mesh)
    else:
        raise NotImplementedError(f"Method {cfg.method} not implemented for DOS calculation.")
    fly.save_data(outdir + prefix + '_DOS.h5', dos, [], [[]], w_pts)
    print("Saved to ", outdir + prefix + "_DOS.h5")
