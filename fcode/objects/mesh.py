import fcode.config as cfg
import numpy as np

def energy_mesh():
    BZ = np.array(cfg.brillouin_zone)
    k1, k2, k3 = np.meshgrid(np.arange(nk1)/nk1, np.arange(nk2)/nk2, np.arange(nk3)/nk3)
    k_points = np.stack((k1.ravel(), k2.ravel(), k3.ravel()), axis=-1)
    k_points = k_points @ BZ.T
    e_pts = fcode.epsilon(1, k_points)
    e_mesh = e_pts.reshape(nk1, nk2, nk3)
    return e_mesh

