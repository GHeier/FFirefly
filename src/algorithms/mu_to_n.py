import firefly
import firefly.config as cfg

import numpy as np

outdir = cfg.outdir
prefix = cfg.prefix
dim = cfg.dimension
nbnd = cfg.nbnd
BZ = np.array(cfg.brillouin_zone)

def get_bandwidth():
    nk1, nk2, nk3 = [200, 200, 200]
    if (dim == 2):
        nk3 = 1
    k1, k2, k3 = np.meshgrid(np.arange(nk1)/nk1, np.arange(nk2)/nk2, np.arange(nk3)/nk3)
    k_points = np.stack((k1.ravel(), k2.ravel(), k3.ravel()), axis=-1)
    k_points = k_points @ BZ.T

    band = firefly.Bands()

    min_e, max_e = [1000, -1000]
    for n in range(1, nbnd):
        e_pts = band(n, k_points)
        e_mesh = e_pts.reshape(nk1, nk2, nk3)
        min_e = min(np.min(e_mesh), min_e)
        max_e = max(np.max(e_mesh), max_e)
    return float(min_e), float(max_e)


def mu_vs_n():
    dos = firefly.Field_R(outdir + prefix + "_DOS.dat")
    min_e, max_e = get_bandwidth()
    bandwidth = max_e - min_e
    de = bandwidth / 1000
    print("Bandwidth: ", bandwidth)
    num = 0
    e_list = list()
    num_list = list()
    for i in range(1000):
        e = min_e + de * i
        num += 2 * de * dos(e) # 2 is for spin degeneracy
        e_list.append(e)
        num_list.append(num)
    save(e_list, num_list)

def save(x, y):
    data = np.column_stack((y, x))
    np.savetxt(outdir + prefix + "_N_vs_mu.dat", data, header="# w f", comments='', fmt='%.6f')
    print("Saved to ", outdir + prefix + "_N_vs_mu.dat")
        
