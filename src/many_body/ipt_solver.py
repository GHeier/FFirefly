from triqs.gf import *
import numpy as np
from math import pi
import triqs.gf.dlr_crm_dyson_solver
from triqs.dos import DOSFromFunction, HilbertTransform
import firefly.config as cfg

def iw_to_tau_dlr(Giw_dlr):
    G_dlr_coeff = make_gf_dlr(Giw_dlr)
    Gtau_dlr = make_gf_dlr_imtime(G_dlr_coeff)
    return Gtau_dlr

def tau_to_iw_dlr(Gtau_dlr):
    G_dlr = make_gf_dlr(Gtau_dlr)
    G_iw = make_gf_dlr_imfreq(G_dlr)
    return G_iw

def iw_dlr_to_w(G_iw_dlr, beta, w_min=-5.0, w_max=5.0, n_w=500):
    n_iw_standard = 100
    imfreq_mesh = MeshImFreq(beta=beta, statistic='Fermion', n_iw=n_iw_standard)
    Giw_standard = Gf(mesh=imfreq_mesh, target_shape=[])

    # Sample DLR Green's function on standard mesh
    G_dlr_coeff = make_gf_dlr(G_iw_dlr)
    Giw_temp = make_gf_imfreq(G_dlr_coeff, n_iw=n_iw_standard)

    # Create real frequency mesh
    Gw = Gf(mesh=MeshReFreq(window = (w_min, w_max), n_w=n_w), target_shape=[1,1])

    # Perform Pade continuation
    Gw.set_from_pade(Giw_temp)
    return Gw

class IPTSolver:
    def __init__(self, beta, H, n_loops=100, mix=0.10, tol=1e-6, w_max=1.2*4, eps=1e-14):
        self.beta = beta
        self.H = H
        self.max_loops = n_loops
        self.mix = mix
        self.tol = tol
        self.mu = cfg.fermi_energy

        # Matsubara frequency Green's functions
        dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        self.G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
        self.Sigma_iw = self.G_iw.copy()
        self.Sigma_iw.zero()
        self.G_iw << H(Sigma = self.Sigma_iw, mu=self.mu)
        print(self.G_iw)
        self.G0_iw = self.G_iw.copy() # self.G0 will be set by the user after initialization
        print(self.G_iw.data.shape)

        # Imaginary time
        tau_mesh = MeshDLRImTime(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        self.G_tau = Gf(mesh=tau_mesh, target_shape=[1,1])
        self.Sigma_tau = self.G_tau.copy()

    def solve(self, U):
        self.G_tau = iw_to_tau_dlr(self.G0_iw)
        self.Sigma_tau << (U**2) * self.G_tau * self.G_tau * self.G_tau
        Sigma_iw = tau_to_iw_dlr(self.Sigma_tau)
        self.Sigma_iw = self.mix * Sigma_iw + (1.0 - self.mix) * self.Sigma_iw

        # Dyson
        #self.G_iw << inverse(inverse(self.G0_iw) - self.Sigma_iw)
        self.G_iw << self.H(Sigma=self.Sigma_iw, mu=self.mu)
        #self.G0_iw << inverse( iOmega_n - t**2 * self.G_iw )
        self.G0_iw << inverse(inverse(self.G_iw) + self.Sigma_iw)
        #self.G_iw = self.G0_iw * self.mix + self.G_iw * (1.0 - self.mix)

    def solve_bethe_lattice(self, U):
        self.G_tau = iw_to_tau_dlr(self.G0_iw)
        self.Sigma_tau << (U**2) * self.G_tau * self.G_tau * self.G_tau
        self.Sigma_iw = tau_to_iw_dlr(self.Sigma_tau)

        self.G_iw = inverse(inverse(self.G0_iw) - self.Sigma_iw)
        self.G0_iw << inverse( iOmega_n - t**2 * self.G_iw )
        #self.G_iw = self.G0_iw * self.mix + self.G_iw * (1.0 - self.mix)

    def loop(self, U, bethe_lattice=False):
        for i in range(self.max_loops):
            G_iw_old = self.G_iw.copy()
            if bethe_lattice:
                self.solve_bethe_lattice(U)
            else:
                self.solve(U)
            err = abs(self.G_iw.data - G_iw_old.data).max()
            print("IPT loop %d, err = %.3e" % (i+1, err))
            if err < self.tol:
                break


#from triqs.plot.mpl_interface import *
## change scale of all figures to make them bigger
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#
#t = 1.0
#U = 4.0
#beta = 50
#
## Define Bethe lattice semicircular DOS: rho(e) = sqrt(4*t^2 - e^2) / (2*pi*t^2)
#def bethe_dos(e):
#    """Semicircular density of states for Bethe lattice with half-bandwidth 2*t"""
#    half_bandwidth = 2*t
#    if abs(e) < half_bandwidth:
#        return np.sqrt(half_bandwidth**2 - e**2) / (np.pi * half_bandwidth)
#    else:
#        return 0.0
#
## Create DOS object from function
#dos = DOSFromFunction(function=bethe_dos, x_min=-2*t, x_max=2*t, n_pts=1000, name='bethe_dos')
#
## Create Hilbert transform object
#H = HilbertTransform(dos)
#
#
#def test():
#    S = IPTSolver(beta = beta, H = H)
#    S.G_iw << SemiCircular(2*t)
#
#    G0_tau = iw_to_tau_dlr(S.G_iw)
#    Sigma_tau = G0_tau.copy()
#    Sigma_tau << (U**2) * G0_tau * G0_tau * G0_tau
#    Sigma_iw = tau_to_iw_dlr(Sigma_tau)
#
#    # Dyson
#    G_iw = inverse(inverse(S.G_iw) - Sigma_iw)
#
#    S.G0_iw << SemiCircular(2*t)
#    S.solve_bethe_lattice(U = U)
#    bethe_G = S.G_iw.data[:, 0, 0].imag
#    bethe_G_w = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#    S.G0_iw << SemiCircular(2*t)
#    S.G_iw << SemiCircular(2*t)
#    S.solve(U = U)
#    gen_G = S.G_iw.data[:, 0, 0].imag
#    gen_G_w = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#
#    plt.plot(gen_G)
#    plt.plot(bethe_G)
#    plt.show()
#
#    oplot(-gen_G_w.imag/pi, label='General Lattice')
#    oplot(-bethe_G_w.imag/pi, label='Bethe Lattice')
#    plt.show()
#
#
#def loop_comparison():
#    S = IPTSolver(beta = beta, H = H)
#    S.G0_iw << SemiCircular(2*t)
#    S.G_iw << SemiCircular(2*t)
#    S.loop(U=U)
#    Gw = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#    oplot(-Gw.imag/pi, label='General Lattice')
#
#    S.G0_iw << SemiCircular(2*t)
#    S.G_iw << SemiCircular(2*t)
#    S.loop(U=U, bethe_lattice=True)
#    Gw = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#    oplot(-Gw.imag/pi, label='Bethe Lattice')
#
#    plt.ylim(0,0.35)
#    plt.show()
#
#def mit_test():
#    fig = plt.figure(figsize=(12,8))
#
#    S = IPTSolver(beta = beta, H = H)
#    S.G0_iw << SemiCircular(2*t)
#    S.G_iw << SemiCircular(2*t)
#    S.loop(U=U)
#    Gw = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#    oplot(-Gw.imag/pi, label='Bethe Lattice')
#    S.G0_iw << SemiCircular(2*t)
#    S.G_iw << SemiCircular(2*t)
#    for i in range(S.max_loops):
#        S.solve(U = U)
#
#        Gw = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#
#        if i % 5 == 0:
#            oplot(-Gw.imag/pi, figure = fig, label = "Iteration = %i" % (i+1), name=r"$\rho$")
#
#    #oplot(-Gw.imag/pi, label='Bethe Lattice')
#    plt.ylim(0,0.35)
#    plt.show()
#
#def mit():
#    fig = plt.figure(figsize=(12,8))
#    U_list = [0, 2, 4, 6, 8]
#    #U_list = [0]
#    n_U = len(U_list)
#    pn = 0 # iteration counter for plotting
#    S = IPTSolver(beta = beta, H = H)
#    for U in U_list:
#        S.G0_iw << SemiCircular(2*t)
#        S.G_iw << SemiCircular(2*t)
#        S.loop(U=U)
#        # Get the real-axis with Pade approximation
#        G_w = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#
#        # plotting
#        ax = fig.add_axes([0,1.-(pn+1)/n_U,1,1./n_U]) # subplot
#        ax.set_xticklabels([])
#        ax.set_yticklabels([])
#        oplot(-G_w.imag/pi, linewidth=1, label = "U = %.2f" % U)
#        plt.xlim(-8,8)
#        plt.ylim(0,0.35)
#        plt.ylabel("")
#        pn = pn + 1
#    plt.show()
#
#def plot_G_iw(G_iw):
#    iwn = G_iw.copy()
#    iwn << iOmega_n
#    start = int(len(iwn.data) / 2)
#    #print(start)
#    plt.plot(iwn.data.imag[start:,0,0], G_iw.data.imag[start:,0,0])
#    plt.xlabel('iwn')
#    plt.ylabel('G(iwn)')
#    plt.show()

#test()
#loop_comparison()
#mit_test()
#mit()

#S = IPTSolver(beta = beta, H = H)
#S.G0_iw << SemiCircular(2*t)
#S.loop(U=3)
#plot_G_iw(S.G_iw)

#S.G0_iw << iOmega_n
#
#size = len(S.G0_iw.data)
#for i in range(size):
#    n = i - size / 2
#    print((2*n+1)*pi / beta, S.G0_iw.data[i][0][0].imag)

#fig = plt.figure(figsize=(12,8))

#S = IPTSolver(beta = beta, H = H)
#S.G0_iw << SemiCircular(2*t)
#S.G_iw << SemiCircular(2*t)
#S.loop(U=U, bethe_lattice=False)
#Gw = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#oplot(-Gw.imag/pi, label='Bethe Lattice')
#plt.show()

#fig = plt.figure(figsize=(12,8))
#S = IPTSolver(beta = beta, H = H)
#S.G0_iw << SemiCircular(2*t)
#S.G_iw << SemiCircular(2*t)
#
#for i in range(S.max_loops):
#    S.solve(U = U)
#    #S.G0_iw = S.G_iw.copy()
#
#    Gw = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#
#    if i % 8 == 0:
#        oplot(-Gw.imag/pi, figure = fig, label = "Iteration = %i" % (i+1), name=r"$\rho$")
#plt.ylim(0,0.35)
#plt.show()


#fig = plt.figure(figsize=(6,6))
#pn = 0 # iteration counter for plotting
#
#for U in [0, 2, 3, 4, 5, 6, 7]:
#
#    S = IPTSolver(beta = beta)
#    S.G_iw << SemiCircular(2*t)
#
#    # DMFT
#    for i in range(n_loops):
#        S.G0_iw << inverse( iOmega_n - t**2 * S.G_iw )
#        S.solve(U)
#
#    # Get the real-axis with Pade approximation
#    G_w = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#
#    # plotting
#    ax = fig.add_axes([0,1.-(pn+1)/6.,1,1./6.]) # subplot
#    ax.set_xticklabels([])
#    ax.set_yticklabels([])
#    oplot(-G_w.imag/pi, linewidth=1, label = "U = %.2f" % U)
#    plt.xlim(-8,8)
#    plt.ylim(0,0.35)
#    plt.ylabel("")
#    pn = pn + 1
#plt.show()

