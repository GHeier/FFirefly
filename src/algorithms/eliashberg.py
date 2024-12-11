from ..config.load.python_config import *
from fcode import *

class Eliashberg:
    def __init__(self, k_pts, scf_tol=1e-4, maxiter=100, mix=0.2, verbose=True):
        self.k_pts = k_pts
        self.scf_tol = scf_tol
        self.maxiter = maxiter
        self.mix = mix
        self.verbose = verbose

        beta = 1 / Temperature
        IR_tol = 1e-10
        wmax = 2*get_bandwidth(nk1, nk2, nk3)
        IR_basis_set = sparse_ir.FiniteTempBasisSet(beta, wmax, eps=IR_tol)
        self.phi = np.ones((len(k_pts), len(k_pts), 1, pts_matsubara), dtype=complex)
        self.z = [1] * len(k_pts)
        self.chi = [1] * len(k_pts)

    def get_denominator(k1, w1, phi_el, Z_el, chi_el):
        term1 = w1 * Z_el
        term2 = epsilon(k1) + chi_el
        term3 = phi_el

        term1 = term1.real**2 + term1.imag**2
        term2 = term2.real**2 + term2.imag**2
        term3 = term3.real**2 + term3.imag**2

        val = term1 + term2 + term3
        return val

    def k_integral(k, w, w1, phi, Z, chi, arr, dim1, dim2, dim3, dim4, xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax):
        phi_int = 0
        Z_int = 0
        chi_int = 0
        n = int(np.round((w.imag / (np.pi / beta) - 1.0) / 2.0))
        m_z = pts_integrate
        if dim == 2:
            m_z = 1

        for i in range(1, pts_k + 1):
            for j in range(1, pts_k + 1):
                for l in range(1, m_z + 1):
                    k1 = Vec(-np.pi + i * 2 * np.pi / pts_k, -np.pi + j * 2 * np.pi / pts_k, -np.pi + l * 2 * np.pi / pts_k)

                    phi_el = phi[i-1, j-1, l-1, n]
                    Z_el = Z[i-1, j-1, l-1, n]
                    chi_el = chi[i-1, j-1, l-1, n]

                    V1 = FLEX_interp(k - k1, w - w1, arr, xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax)
                    V2 = FLEX_interp(k + k1, w - w1, arr, xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax)
                    V_even = 0.5 * (V1 + V2)

                    denom = get_denominator(k1, w1, phi_el, Z_el, chi_el)
                    phi_int += V_even * phi_el / denom
                    Z_int += V_even * w1 / w * Z_el / denom
                    chi_int += V_even * (epsilon(k1) + chi_el) / denom

        phi_int /= ((2 * np.pi) ** dim * pts_k ** 2 * m_z)
        Z_int /= ((2 * np.pi) ** dim * pts_k ** 2 * m_z)
        chi_int /= ((2 * np.pi) ** dim * pts_k ** 2 * m_z)

    def eliashberg_sum(self):
        self.phi = 
        for a in range(1, pts_k + 1):
            print(a, "out of ", pts_k)
            for b in range(1, pts_k + 1):
                for c in range(1, m_z + 1):
                    k = Vec(-np.pi + a * 2 * np.pi / pts_k, -np.pi + b * 2 * np.pi / pts_k, -np.pi + c * 2 * np.pi / pts_k)

                    for i in range(1, pts_matsubara + 1):
                        w = complex(0, (2 * i + 1) * np.pi / beta)
                        for j in range(1, pts_matsubara + 1):
                            w1 = complex(0, (2 * j + 1) * np.pi / beta + 0.0001)
                            k_integral(k, w, w1, phi, Z, chi, arr, dim1, dim2, dim3, dim4, xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax)
                            new_phi[a-1, b-1, c-1, i-1] += phi_int
                            new_Z[a-1, b-1, c-1, i-1] += Z_int
                            new_chi[a-1, b-1, c-1, i-1] += chi_int

                        new_phi[a-1, b-1, c-1, i-1] /= beta
                        new_Z[a-1, b-1, c-1, i-1] = 1 + 1 / beta * new_Z[a-1, b-1, c-1, i-1]
                        new_chi[a-1, b-1, c-1, i-1] /= -beta

    def evaluate_eliashberg():
        phi = np.ones((pts_k, pts_k, 1, pts_matsubara), dtype=complex)
        Z = np.ones((pts_k, pts_k, 1, pts_matsubara), dtype=complex)
        chi = np.ones((pts_k, pts_k, 1, pts_matsubara), dtype=complex)
        new_phi = np.zeros_like(phi)
        new_Z = np.zeros_like(Z)
        new_chi = np.zeros_like(chi)

        arr = read_4D_data("/home/g/Research/bcs_diagonalization/matsubara_cube.dat")
        dim1, dim2, dim3, dim4 = arr.shape
        xmin, xmax = 0.0, 2 * np.pi
        ymin, ymax = 0.0, 2 * np.pi
        zmin, zmax = 0.0, 2 * np.pi
        wmin, wmax = 0.0, 10.45

        print("Starting Eliashberg calculation")
        for i in range(10):
            eliashberg_sum(phi, Z, chi, arr, dim1, dim2, dim3, dim4, xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax)
            phi = new_phi
            Z = new_Z
            chi = new_chi
            print("(Sample) Phi: ", phi[0, 0, 0, 0], " Z: ", Z[0, 0, 0, 0], " Chi: ", chi[0, 0, 0, 0])
            print("(Max) Phi: ", np.max(np.real(phi)), " Z: ", np.max(np.real(Z)), " Chi: ", np.max(np.real(chi)))
            filename = f"phi_{i}.dat"
            save_4D_data(phi, filename)
            filename = f"Z_{i}.dat"
            save_4D_data(Z, filename)
            filename = f"chi_{i}.dat"
            save_4D_data(chi, filename)

