from triqs.lattice import BravaisLattice, BrillouinZone
from triqs.gf import Gf, MeshProduct, MeshBrZone, MeshImFreq, MeshReFreq
import numpy as np
from math import cos, pi

from triqs.plot.mpl_interface import oplot,plt
from triqs.gf import *
#from triqs.gf import Gf, SemiCircular
#from triqs.gf import MeshReFreq, MeshImFreq
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
#from triqs.gf import make_gf_dlr, make_gf_dlr_imtime, make_gf_imfreq, make_gf_imtime

beta = 50
mu = 0.0
t = 1
U = 4
D = 4*t # Bandwidth

BL = BravaisLattice([(1,0,0), (0,1,0)]) # Two unit vectors in R3
BZ = BrillouinZone(BL)

# n_k denotes the number of k-points for each dimension
n_k = 128
k_mesh = MeshBrZone(bz=BZ, n_k=n_k)

iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=1.2*D, eps=1e-14)
k_iw_mesh = MeshProduct(k_mesh, iw_mesh)

# Recall that for an empty target_shape G0 has values that are scalars instead of matrices.
Gw = Gf(mesh=k_iw_mesh, target_shape=[])

def eps(k):
    return -2*t * (cos(k[0]) + cos(k[1]))

iw_arr = np.array(list(iw_mesh.values()))
k_arr  = np.array(list(k_mesh.values()))
np_eps = np.vectorize(eps, signature='(d)->()')

eps_arr = np_eps(k_arr)
Gw.data[:] = 1.0 / (iw_arr[None,::] + mu - eps_arr[::,None])

Gw0 = Gw.copy()
Ew = Gw.copy()
tau_mesh = MeshDLRImTime(beta=beta, statistic='Fermion', w_max = 1.2*D, eps=1e-14)
Gt = Gf(mesh=tau_mesh, target_shape=[1,1])
Et = Gt.copy()

def IPT_iteration(Gw0, Gw, Gt, Ew, Et, U):
    Gt << Fourier(Gw0)
    Et << (U**2) * Gt * Gt * Gt
    Ew << Fourier(Et)

    # Dyson
    Gw << inverse(inverse(Gw0) - Ew)


def loop(Gw, Gt, Ew, Et, U):
    iters = 0
    err = float('inf')
    while err > 1e-4 and iters < 100:
        IPT_iteration(Gw0, Gw, Gt, Ew, Et, U)
        err = np.linalg.norm((G1 - G2).data)
        iters += 1

loop(Gw, Gt, Ew, Et, U)

#-----------------------------------------------------------------------------------------------
#dlr_iw_mesh = MeshDLRImFreq(beta=50, statistic='Fermion', w_max=1.2, eps=1e-14)
#Giw_dlr = Gf(mesh= dlr_iw_mesh, target_shape=[1,1])
#Giw_dlr << SemiCircular(1.0)
#
#iw_mesh = MeshImFreq(beta=50, S='Fermion', n_iw=100)
#Giw = Gf(mesh=iw_mesh, target_shape=[1,1])
#Giw << SemiCircular(1.0)
#
#w_mesh = MeshReFreq(window=(-4,4), n_w=500)
#Gw = Gf(mesh=w_mesh, target_shape=[1,1])
#Giw_from_dlr = make_gf_imfreq(Giw_dlr, n_iw=100)
#Gw.set_from_pade(Giw_from_dlr)
#
#oplot(-Gw.imag/pi, name=r"$\rho$")
#plt.show()
