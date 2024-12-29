# Get fmodule
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))
from ..fmodule import *
from ..config.load import config as cfg

import numpy as np

def energy_mesh(pts):
    nk1, nk2, nk3 = pts
    BZ = np.array(cfg.brillouin_zone)
    k1, k2, k3 = np.meshgrid(np.arange(nk1)/nk1, np.arange(nk2)/nk2, np.arange(nk3)/nk3)
    k_points = np.stack((k1.ravel(), k2.ravel(), k3.ravel()), axis=-1)
    k_points = k_points @ BZ.T
    e_pts = epsilon(1, k_points)
    e_mesh = e_pts.reshape(nk1, nk2, nk3)
    return e_mesh

