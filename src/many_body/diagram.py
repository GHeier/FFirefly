from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
from triqs.gf import MeshProduct, MeshBrillouinZone
from triqs_tprf.lattice import fourier_tr_to_wr, fourier_wk_to_wr, fourier_wr_to_tr, fourier_wr_to_wk, chi_wr_from_chi_tr, chi_wk_from_chi_wr, chi_tr_from_chi_wr, chi_wr_from_chi_wk
import numpy as np

import firefly as fly
import firefly.config as cfg

class Diagram:
    def __init__(self, obj, statistic):
        varspace = describe_mesh(obj)
        if varspace not in ['wk', 'tr', 'w', 't']:
            print(f"Diagram initialized in {varspace} space.")
            raise ValueError("Object must be MeshDLRImFreq or MeshDLRImTime")

        self.statistic = statistic
        self.shape = obj.data.shape
        self.nw = self.shape[0]
        if varspace != 'w' and varspace != 't':
            self.nk = self.shape[1]
        self.ind_dim = len(self.shape) - 2

        if varspace == 'wk':
            self.obj_wk = obj
            self.wk_to_tr()
        elif varspace == 'tr':
            self.obj_tr = obj
            self.tr_to_wk()
        elif varspace == 'w':
            self.obj_w = obj
            self.w_to_t()
        elif varspce == 't':
            self.obj_t = obj
            self.t_to_w()
        else:
            raise ValueError("varspace must be 'wk' or 'tr' or 'w' or 't'")

        # Extract w-points from mesh (Matsubara frequencies)
        # Must be done AFTER transformation to wk space
        mesh_w = obj.mesh.components[0]
        # For Matsubara frequencies, use imaginary part
        self.w_points = np.array([float(iw.imag) for iw in mesh_w], dtype=np.float32)

    def wk_to_tr(self):
        if self.statistic == 'Fermion':
            self.obj_wr = fourier_wk_to_wr(self.obj_wk)
            self.obj_tr = fourier_wr_to_tr(self.obj_wr)
        else:
            self.obj_wr = chi_wr_from_chi_wk(self.obj_wk)
            self.obj_tr = chi_tr_from_chi_wr(self.obj_wr)

    def tr_to_wk(self):
        if self.statistic == 'Fermion':
            self.obj_wr = fourier_tr_to_wr(self.obj_tr)
            self.obj_wk = fourier_wr_to_wk(self.obj_wr)
        else:
            self.obj_wr = chi_wr_from_chi_tr(self.obj_tr, nw=self.nw)
            self.obj_wk = chi_wk_from_chi_wr(self.obj_wr)

    def w_to_t(self):
        self.dlr = make_gf_dlr(self.obj_w)
        self.obj_t = make_gf_dlr_imtime(self.dlr)

    def tau_to_iw_dlr(self):
        self.dlr = make_gf_dlr(self.obj_t)
        self.obj_w = make_gf_dlr_imfreq(self.dlr)

    def save(self, filename):
        mesh, BZ = extract_mesh_and_bz(self.obj_wk)
        # Reshape data to match mesh dimensions (nw, nkx, nky, nkz)
        obj = np.reshape(self.obj_wk.data, mesh)
        fly.save_data(filename, obj.T, mesh=np.array(mesh, dtype=np.int32), domain=BZ, w_points=self.w_points)

    def save_as_w(self, filename):
        obj_w = np.sum(self.obj_wk.data, axis=1) / self.nk  # Sum over k-points
        obj_w = np.reshape(obj_w, (self.nw, ))
        fly.save_data(filename, obj_w, mesh=None, domain=None, w_points=self.w_points)


def copy(diagram):
    new_diag = Diagram(diagram.obj_wk, diagram.statistic)
    return new_diag

def dot(diagram1, diagram2):
    # Setup for one-band only for now
    new_obj = copy(diagram1)
    #print(diagram1.obj_tr.data.shape)
    #print(diagram2.obj_tr.data[:, :, 0, 0].shape)
    new_obj.obj_tr.data[:] = diagram1.obj_tr.data * diagram2.obj_tr.data[:, :, 0, 0]
    return new_obj

def describe_mesh(G):
    mesh = G.mesh
    if isinstance(mesh, MeshProduct):
        m0, m1 = mesh.components

        if isinstance(m0, MeshDLRImFreq):
            return "wk"
        elif isinstance(m0, MeshDLRImTime):
            return "tr"
        else:
            return "unknown MeshProduct components"

    elif isinstance(mesh, MeshDLRImFreq):
        return "w"
    elif isinstance(mesh, MeshDLRImTime):
        return "t"
    elif isinstance(mesh, BrillouinZoneMesh):
        return "k"
    else:
        return "unknown mesh type"


def extract_mesh_and_bz(G):
    """
    Extract mesh dimensions and Brillouin zone from a Green's function.

    Args:
        G: TRIQS Green's function object with MeshProduct(MeshDLRImFreq/MeshDLRImTime, MeshBrillouinZone)

    Returns:
        mesh: tuple of mesh dimensions (nw, nk_x, nk_y, nk_z) or (nw, nk_x, nk_y) for 2D
        brillouin_zone: numpy array of shape (dim, dim) containing reciprocal lattice vectors

    Example:
        For a 2D square lattice with 60x60 k-mesh:
        mesh = (28, 60, 60, 1)
        brillouin_zone = [[6.28, 0.0, 0.0], [0.0, 6.28, 0.0], [0.0, 0.0, 6.28]]
    """
    mesh = G.mesh

    if not isinstance(mesh, MeshProduct):
        raise ValueError("Green's function must have a MeshProduct mesh")

    mesh_w, mesh_k = mesh.components

    # Get frequency mesh size
    nw = len(mesh_w)

    # Get k-mesh dimensions
    if not isinstance(mesh_k, MeshBrillouinZone):
        raise ValueError("Second mesh component must be MeshBrillouinZone")

    # Get k-mesh size - dims attribute gives (nkx, nky, nkz)
    k_dims = mesh_k.dims

    # Combine into full mesh tuple (nw, nkx, nky, nkz)
    mesh_tuple = (nw,) + tuple(k_dims)

    # Extract Brillouin zone (reciprocal lattice vectors)
    # The bz.units attribute contains the reciprocal lattice vectors
    bz_matrix = mesh_k.bz.units

    # Convert to numpy array (shape: 3x3)
    brillouin_zone = np.array(bz_matrix, dtype=np.float32)

    return mesh_tuple, brillouin_zone

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

