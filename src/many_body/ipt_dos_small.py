from triqs.gf import *
import numpy as np
from math import pi
import triqs.gf.dlr_crm_dyson_solver
from triqs.dos import DOSFromFunction, HilbertTransform

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

def p(e):
    """Bethe lattice semicircular density of states for half-bandwidth 2"""
    if abs(e) < 2.0:
        return np.sqrt(4.0 - e**2) / (2.0 * pi)
    else:
        return 0.0

def integrate_de(Sigma, mu=0.0):
    """Integrate over energy using Bethe lattice semicircular DOS"""
    n_pts = 1000
    e_vals = np.linspace(-2.0, 2.0, n_pts)
    de = e_vals[1] - e_vals[0]
    pe_vals = np.vectorize(p)(e_vals)
    integral = 0.0
    iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
    iw << iOmega_n
    for e in e_vals:
        integral += p(e) / (iw - e - Sigma)
    integral *= de
    return integral

class IPTSolver:
    def __init__(self, beta, w_max=1.2*4, eps=1e-14):
        self.beta = beta

        # Matsubara frequency Green's functions
        dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        self.G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
        self.G0_iw = self.G_iw.copy() # self.G0 will be set by the user after initialization
        self.Sigma_iw = self.G_iw.copy()

        # Imaginary time
        tau_mesh = MeshDLRImTime(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        self.G0_tau = Gf(mesh=tau_mesh, target_shape=[1,1])
        self.Sigma_tau = self.G0_tau.copy()

    def solve(self, U, hilbert_transform=None):
        self.G0_tau = iw_to_tau_dlr(self.G0_iw)
        self.Sigma_tau << (U**2) * self.G0_tau * self.G0_tau * self.G0_tau
        self.Sigma_iw = tau_to_iw_dlr(self.Sigma_tau)

        # Use Hilbert transform if provided, otherwise use Dyson equation
        if hilbert_transform is not None:
            # Compute lattice Green's function using full DOS integral
            # G(iw) = ∫ rho(e) de / (iw - e - Sigma(iw))
            G_iw = hilbert_transform(Sigma=self.Sigma_iw, mu=0.0)
        else:
            # Dyson equation
            G_iw = inverse(inverse(self.G0_iw) - self.Sigma_iw)

        self.G_iw = G_iw * mix + self.G_iw * (1.0 - mix)

from triqs.plot.mpl_interface import *
# change scale of all figures to make them bigger
import matplotlib as mpl
import matplotlib.pyplot as plt

t = 1.0
U = 4.0
beta = 50
n_loops = 100
mix = 0.05

# Define Bethe lattice semicircular DOS: rho(e) = sqrt(4*t^2 - e^2) / (2*pi*t^2)
def bethe_dos(e):
    """Semicircular density of states for Bethe lattice with half-bandwidth 2*t"""
    half_bandwidth = 2*t
    if abs(e) < half_bandwidth:
        return np.sqrt(half_bandwidth**2 - e**2) / (np.pi * half_bandwidth)
    else:
        return 0.0

# Create DOS object from function
dos = DOSFromFunction(function=bethe_dos, x_min=-2*t, x_max=2*t, n_pts=1000, name='bethe_dos')

# Create Hilbert transform object
H = HilbertTransform(dos)

S = IPTSolver(beta = beta)
S.G_iw << SemiCircular(2*t)
#S.G0_iw << iOmega_n
#
#size = len(S.G0_iw.data)
#for i in range(size):
#    n = i - size / 2
#    print((2*n+1)*pi / beta, S.G0_iw.data[i][0][0].imag)

#fig = plt.figure(figsize=(12,8))
#
#for i in range(n_loops):
#    # Updated: Use full DOS integral instead of simplified Bethe lattice formula
#    S.G0_iw << inverse(inverse(S.G_iw) + S.Sigma_iw)
#    #S.G0_iw = S.G_iw.copy()
#    S.solve(U=U, hilbert_transform=H)
#
#    Gw = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)
#
#    if i % 20 == 0:
#        oplot(-Gw.imag/pi, figure = fig, label = "Iteration = %i" % (i+1), name=r"$\rho$")
#plt.ylim(0,0.35)
#plt.show()


fig = plt.figure(figsize=(6,6))
pn = 0 # iteration counter for plotting

for U in [0, 2, 3, 4, 5, 6, 7]:

    S = IPTSolver(beta = beta)
    S.G_iw << SemiCircular(2*t)

    # DMFT
    for i in range(n_loops):
        G_old = S.G_iw.copy()
        # Self-consistency: G0^{-1} = G^{-1} + Sigma
        # This replaces the simplified Bethe lattice formula: G0^{-1} = iw - t^2 * G
        S.G0_iw << inverse(inverse(S.G_iw) + S.Sigma_iw)
        S.solve(U, hilbert_transform=H)
        err = np.linalg.norm((G_old - S.G_iw).data)
        if err < 1e-4:
            print("Converged for U = %.2f after %i iterations" % (U, i+1))
            break

    # Get the real-axis with Pade approximation
    G_w = iw_dlr_to_w(S.G_iw, beta, w_min=-8.0, w_max=8.0, n_w=1000)

    # plotting
    ax = fig.add_axes([0,1.-(pn+1)/6.,1,1./6.]) # subplot
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    oplot(-G_w.imag/pi, linewidth=1, label = "U = %.2f" % U)
    plt.xlim(-8,8)
    plt.ylim(0,0.35)
    plt.ylabel("")
    pn = pn + 1
plt.show()

