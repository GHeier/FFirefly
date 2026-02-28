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
smearing = cfg.smearing

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

def get_electron_number(dos, w_pts):
    print("Calculating electron number vs energy...")
    n_list = np.zeros(len(w_pts))
    dw = w_pts[1] - w_pts[0]
    for i in range(len(w_pts)):
        n_list[i] = 2*np.sum(dos[:i]) * dw
    fly.save_data(outdir + prefix + '_E_vs_n.h5', w_pts, [], [[]], n_list)
    print("Saved to ", outdir + prefix + "_E_vs_n.h5")


def run():
    print("Calculating DOS using Gaussian smearing...")
    print(f"Smearing width: {smearing}")
    print(f"Memory estimate: ~{ne*nx*ny*nz*8/1e9:.2f} GB for energy mesh")
    e_mesh = get_e_mesh()
    emin = np.min(e_mesh) + 2 * smearing
    emax = np.max(e_mesh) - 2 * smearing
    grid = np.linspace(emin, emax, ne)
    print(f"Energy range: {emin:.4f} to {emax:.4f} with {ne} points at a spacing of {(emax-emin)/(ne-1):.4f}")
    X = grid[:, None] - e_mesh[None, :]
    dos = np.sum(np.exp(-X**2 / (2 * smearing**2)), axis=1)

    area = np.trapezoid(dos, grid)
    if area > 0:
        dos /= area
    else:
        print("DOS area is zero or negative.")
    print("DOS calculation completed. Saving")
    fly.save_data(outdir + prefix + '_DOS.h5', dos, [], [[]], grid)
    print("Saved to ", outdir + prefix + "_DOS.h5")
    get_electron_number(dos, grid)

if __name__ == "__main__":
    run()
