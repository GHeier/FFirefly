import fcode
import fcode.config as cfg

import numpy as np
import sparse_ir

BZ = cfg.brillouin_zone
#print("BZ:", BZ)

T = cfg.Temperature
beta = 1 / T
nk1, nk2, nk3 = cfg.k_mesh
#kvec = fcode.Vec(4, 2, 3)
#print("kvec:", kvec.x, kvec.y, kvec.z)

#e = fcode.e_test(kvec)
#print(e)

#print("nbnd:", cfg.nbnd)
#print("Temperature:", cfg.Temperature) 

#e = fcode.epsilon(0, kvec)
#print("e:", e)

t = 1
def epsilon(k):
    return -2*t*(np.cos(k.x) + np.cos(k.y) + np.cos(k.z))

def get_bandwidth():
    k1, k2, k3 = np.meshgrid(np.arange(nk1)/nk1, np.arange(nk2)/nk2, np.arange(nk3)/nk3)
    emin, emax = 1000, -1000
    for i in range(1, 2*nk1):
        for j in range(1, 2*nk2):
            for k in range(1, 2*nk3):
                ivec = [i/nk1, j/nk2, k/nk3]
                kvec = np.dot(ivec, BZ)
                kvec = fcode.Vec(kvec[0], kvec[1], kvec[2])
                ek = epsilon(kvec)
                #print("kvec:", kvec.x, kvec.y, kvec.z)
                #print("ek:", ek)
                emin = min(emin, ek)
                emax = max(emax, ek)
    return emax - emin


def eliashberg():
    W = get_bandwidth()
    W = round(W, 6)
    print("Bandwidth:", W)
    wmax = 2*W
    IR_tol = 1e-5
    IR_basis_set = sparse_ir.FiniteTempBasisSet(beta, wmax, eps=IR_tol)
    print("Loading C config")
    fcode.load_c_config()
    print("Loading C++ config")
    fcode.load_cpp_config()
    print("Testing vec and epsilon")
    test_vec = fcode.Vec(1, 2, 3)
    print("test_vec:", test_vec.x, test_vec.y, test_vec.z)
    a = fcode.get_nbnd()
    print("nbnd:", a)
    fcode.set_nbnd(10)
    a = fcode.get_nbnd()
    print("nbnd:", a)
    b = fcode.get_onsite_U()
    print("onsite_U:", b)
    fcode.set_onsite_U(2)
    print("onsite_U:", b)
    b = fcode.get_onsite_U()
    print("onsite_U:", b)
    fcode.set_FS_only(True)
    d = fcode.get_k_mesh()
    print("k_mesh:", d)
    fcode.set_k_mesh([4, 4, 4])
    print("k_mesh:", d)


    test_epsilon = fcode.epsilon(0, test_vec)
    print("test_epsilon:", test_epsilon)


class Mesh:
    """
    Holding class for k-mesh and sparsely sampled imaginary time 'tau' / Matsubara frequency 'iwn' grids.
    Additionally it defines the Fourier transform routines 'r <-> k'  and 'tau <-> l <-> wn'.
    """
    def __init__(self,IR_basis_set,nk1,nk2,nk3):
        self.IR_basis_set = IR_basis_set

        # generate k-mesh and dispersion
        self.nk1, self.nk2, self.nk3, self.nk = nk1, nk2, nk3, nk1*nk2*nk3
        self.k1, self.k2, self.k3 = np.meshgrid(np.arange(self.nk1)/self.nk1, np.arange(self.nk2)/self.nk2, np.arange(self.nk3)/self.nk3)
        self.ek = -2*t*( np.cos(2*np.pi*self.k1) + np.cos(2*np.pi*self.k2) + np.cos(2*np.pi*self.k3)).reshape(self.nk)

        # lowest Matsubara frequency index
        self.iw0_f = np.where(self.IR_basis_set.wn_f == 1)[0][0]
        self.iw0_b = np.where(self.IR_basis_set.wn_b == 0)[0][0]

        ### Generate a frequency-momentum grid for iwn and ek (in preparation for calculating the Green function)
        # frequency mesh (for Green function)
        self.iwn_f = 1j * self.IR_basis_set.wn_f * np.pi * T
        self.iwn_f_ = np.tensordot(self.iwn_f, np.ones(self.nk), axes=0)

        # ek mesh
        self.ek_ = np.tensordot(np.ones(len(self.iwn_f)), self.ek, axes=0)

    def smpl_obj(self, statistics):
        """ Return sampling object for given statistic """
        smpl_tau = {'F': self.IR_basis_set.smpl_tau_f, 'B': self.IR_basis_set.smpl_tau_b}[statistics]
        smpl_wn  = {'F': self.IR_basis_set.smpl_wn_f,  'B': self.IR_basis_set.smpl_wn_b }[statistics]
        return smpl_tau, smpl_wn

    
    def tau_to_wn(self, statistics, obj_tau):
        """ Fourier transform from tau to in via IR basis """
        smpl_tau, smpl_wn = self.smpl_obj(statistics)

        obj_tau = obj_tau.reshape((smpl_tau.tau.size, self.nk1, self.nk2, self.nk3))
        obj_l   = smpl_tau.fit(obj_tau, axis=0)
        obj_wn  = smpl_wn.evaluate(obj_l, axis=0).reshape((smpl_wn.wn.size, self.nk))
        return obj_wn

    def wn_to_tau(self, statistics, obj_wn):
        """ Fourier transform from iwn to tau via IR basis """
        smpl_tau, smpl_wn = self.smpl_obj(statistics)

        obj_wn  = obj_wn.reshape((smpl_wn.wn.size, self.nk1, self.nk2, self.nk3))
        obj_l   = smpl_wn.fit(obj_wn, axis=0)
        obj_tau = smpl_tau.evaluate(obj_l, axis=0).reshape((smpl_tau.tau.size, self.nk))
        return obj_tau

    
    def k_to_r(self,obj_k):
        """ Fourier transform from k-space to real space """
        obj_k = obj_k.reshape(-1, self.nk1, self.nk2, self.nk3)
        obj_r = np.fft.fftn(obj_k,axes=(1,2))
        obj_r = obj_r.reshape(-1, self.nk)
        return obj_r

    def r_to_k(self,obj_r):
        """ Fourier transform from real space to k-space """
        obj_r = obj_r.reshape(-1, self.nk1, self.nk2, self.nk3)
        obj_k = np.fft.ifftn(obj_r,axes=(1,2))/self.nk
        obj_k = obj_k.reshape(-1, self.nk)
        return obj_k
