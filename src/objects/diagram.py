from triqs.gf import Gf, make_gf_imfreq
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime, MeshImFreq, MeshReFreq
from triqs.gf import MeshProduct, MeshBrillouinZone, make_gf_dlr, make_gf_dlr_imtime, make_gf_dlr_imfreq
from triqs.lattice import BrillouinZone, BravaisLattice
from triqs_tprf.lattice import fourier_tr_to_wr, fourier_wk_to_wr, fourier_wr_to_tr, fourier_wr_to_wk, chi_wr_from_chi_tr, chi_wk_from_chi_wr, chi_tr_from_chi_wr, chi_wr_from_chi_wk
import numpy as np

import firefly as fly
import firefly.config as cfg

class Diagram:
    def __init__(self, obj, statistic):
        """
        Initialize a Diagram from a TRIQS Gf object.

        Args:
            obj: TRIQS Gf object
            statistic: 'Fermion' or 'Boson'
        """
        varspace = describe_mesh(obj)
        if varspace not in ['wk', 'tr', 'w', 't', 'k', 'r']:
            print(f"Diagram initialized in {varspace} space.")
            raise ValueError("Object must be MeshDLRImFreq or MeshDLRImTime")

        self.varspace = varspace
        self.statistic = statistic
        self.shape = obj.data.shape
        self.nw = self.shape[0]
        if varspace == 'wk':
            self.nk = self.shape[1]
        self.ind_dim = len(self.shape) - len(varspace)

        if varspace == 'wk':
            self.obj_wk = obj
            self.wk_to_tr()
        elif varspace == 'tr':
            self.obj_tr = obj
            self.tr_to_wk()
        elif varspace == 'w':
            self.obj_w = obj
            self.w_to_t()
        elif varspace == 't':
            self.obj_t = obj
            self.t_to_w()
        elif varspace == 'k':
            self.nk = self.shape[0]
            self.obj_k = obj
            self.obj_r = obj.copy()
            self.k_to_r()
        elif varspace == 'r':
            self.nk = self.shape[0]
            self.obj_r = obj
            self.obj_k = obj.copy()
            self.r_to_k()
        else:
            raise ValueError("varspace must be 'wk', 'tr', 'w', 't', 'k', or 'r'")

        # Extract w-points from mesh (Matsubara frequencies)
        if varspace == 'wk' or varspace == 'tr':
            # For wk/tr space, mesh is MeshProduct
            mesh_w = self.obj_wk.mesh.components[0]
        else:
            # For w/t space, mesh is directly the frequency/time mesh
            mesh_w = obj.mesh

        # For Matsubara frequencies, use imaginary part
        if isinstance(mesh_w, MeshDLRImFreq):
            self.w_points = np.array([float(iw.imag) for iw in mesh_w], dtype=np.float32)
        else:
            self.w_points = None

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

    def wk_to_wr(self):
        if self.statistic == 'Fermion':
            self.obj_wr = fourier_wk_to_wr(self.obj_wk)
        else:
            self.obj_wr = chi_wr_from_chi_wk(self.obj_wk)

    def wr_to_wk(self):
        if self.statistic == 'Fermion':
            self.obj_wk = fourier_wr_to_wk(self.obj_wr)
        else:
            self.obj_wk = chi_wk_from_chi_wr(self.obj_wr)

    def w_to_t(self):
        self.dlr = make_gf_dlr(self.obj_w)
        self.obj_t = make_gf_dlr_imtime(self.dlr)

    def t_to_w(self):
        self.dlr = make_gf_dlr(self.obj_t)
        self.obj_w = make_gf_dlr_imfreq(self.dlr)

    def k_to_r(self):
        mesh = self.obj_k.mesh.dims
        k_data = np.reshape(self.obj_k.data, mesh)
        r_data = np.fft.ifftn(k_data, axes=tuple(range(1, len(mesh))))
        self.obj_r.data[:] = np.reshape(r_data, self.shape)

    def r_to_k(self):
        mesh, BZ = extract_mesh_and_bz(self.obj_r)
        r_data = np.reshape(self.obj_r.data, mesh)
        k_data = np.fft.fftn(r_data, axes=tuple(range(1, len(mesh))))
        self.obj_k.data[:] = np.reshape(k_data, self.shape)

    def save(self, filename):
        if self.varspace in ['wk', 'tr']:
            mesh, BZ = extract_mesh_and_bz(self.obj_wk)
            # mesh is (nw, nkx, nky, nkz) for k-space dimensions only
            # Append orbital/tensor dimensions for reshaping
            full_shape = mesh + self.shape[2:]  # (nw, nkx, nky, nkz, orbital_dims...)
            obj = np.reshape(self.obj_wk.data, full_shape)
            # Shift k-points to center (only shift spatial dimensions, not orbital indices)
            n_spatial_dims = len(mesh) - 1  # Exclude nw
            obj = np.fft.fftshift(obj, axes=tuple(range(1, n_spatial_dims + 1)))

            # Spatial mesh should ONLY contain k-space dimensions, not tensor indices
            spatial_mesh = np.array(mesh[1:], dtype=np.int32)  # Skip nw: (nkx, nky, nkz)

            n_indices = self.ind_dim  # Number of indices beyond (w, k)
            dim_indices = self.shape[2] if n_indices > 0 else 1  # Size of each index

            # For single-band systems (dim_indices=1), save as scalar field instead of matrix
            if dim_indices == 1 and n_indices > 0:
                n_indices = 0
                dim_indices = 1

            fly.save_data(filename, obj, mesh=spatial_mesh, domain=BZ, w_points=self.w_points,
                         n_indices=n_indices, dim_indices=dim_indices)
        else:
            obj = np.reshape(self.obj_w.data, (self.nw, ))
            fly.save_data(filename, obj, mesh=None, domain=None, w_points=self.w_points)
        print(f"Diagram saved to {filename}")

    def load(self, field):
        data = field.get_data()
        if self.varspace in ['wk', 'tr']:
            mesh, BZ = extract_mesh_and_bz(self.obj_wk)
            mesh = mesh + self.shape[2:]  # Append orbital dimensions
            print("Load Mesh = ", mesh)
            data = np.reshape(data, mesh)  # First reshape to mesh dimensions
            data = np.fft.ifftshift(data, axes=tuple(range(1, len(mesh))))  # Then ifftshift k-axes
            data = np.reshape(data, self.shape)  # Finally reshape to TRIQS data shape
            self.obj_wk.data[:] = data
            self.wk_to_tr()
        else:
            data = np.reshape(data, (self.nw, ))
            self.obj_w.data[:] = data
            self.w_to_t()
        print(f"Diagram loaded from field")

    def save_as_w(self, filename):
        obj_w = np.sum(self.obj_wk.data, axis=1) / self.nk  # Sum over k-points
        obj_w = np.reshape(obj_w, (self.nw, ))
        fly.save_data(filename, obj_w, mesh=None, domain=None, w_points=self.w_points)
        print(f"Diagram(w) saved to {filename}")

    def save_spectral(self, filename):
        beta = self.obj_w.mesh.beta
        wmin, wmax, nw = -5.0, 5.0, 500
        obj_w = iw_dlr_to_w(self.obj_w, beta, w_min=wmin, w_max=wmax, n_w=nw).data
        obj_w = np.reshape(obj_w, (nw, ))
        wpts = np.linspace(wmin, wmax, nw)
        fly.save_data(filename, -obj_w.imag/np.pi, mesh=None, domain=None, w_points=wpts)
        print(f"Diagram(pade) saved to {filename}")

    def init_tail(self):
        spread = 0.01
        tail_vals = spread / (self.w_points ** 2 + spread)
        self.obj_wk.data[:] = tail_vals[:, None, None, None]
        norm = np.sum(self.obj_wk.data * np.conj(self.obj_wk.data)).real
        self.obj_wk.data[:] = self.obj_wk.data / norm
        self.wk_to_tr()

    def copy(self):
        if self.varspace == 'wk':
            data = self.obj_wk.copy()
        elif self.varspace == 'tr':
            data = self.obj_tr.copy()
        elif self.varspace == 't':
            data = self.obj_t.copy()
        elif self.varspace == 'w':
            data = self.obj_w.copy()
        return Diagram(data, self.statistic)

    def zero(self):
        if self.varspace == 'wk' or self.varspace == 'tr':
            self.obj_wk.zero()
            self.obj_tr.zero()
        elif self.varspace == 'w' or self.varspace == 't':
            self.obj_w.zero()
            self.obj_t.zero()

    def fill_from_field(self, field):
        if self.varspace == 'wk' or self.varspace == 'tr':
            BZ = self.obj_wk.mesh.components[1].bz.units
            mesh = self.obj_wk.mesh.components[1].dims
            kpt = np.linspace(0, 1, mesh[0], endpoint=False)
            if len(mesh) == 2:
                kx, ky = np.meshgrid(kpt, kpt, indexing='ij')
                kx = kx.flatten()
                ky = ky.flatten()
                kz = np.zeros_like(kx)
            elif len(mesh) == 3:
                kx, ky, kz = np.meshgrid(kpt, kpt, kpt, indexing='ij')
                kx = kx.flatten()
                ky = ky.flatten()
                kz = kz.flatten()
            else:
                raise ValueError("Mesh must be 2D or 3D")
            kpts = np.vstack((kx, ky, kz)).T @ BZ
            for iw in range(self.nw):
                self.obj_wk.data[iw, :, 0, 0] = field(kpts, self.w_points[iw])
            self.wk_to_tr()
        elif self.varspace == "w":
            for iw in range(self.nw):
                self.obj_w.data[iw] = field(self.w_points[iw])
            self.w_to_t()
        else:
            raise ValueError("fill_diagram_from_field only implemented for 'wk' and 'w' spaces")

def make_local(data):
    return np.einsum('wknm->wnm', data) / data.shape[1]

def get_renorm(loc_sigma, w_points):
    loc_sigma = loc_sigma[:, 0, 0]
    signs = np.sign(loc_sigma.imag)
    diff = np.diff(signs)
    zero_crossings = np.where(diff != 0)[0]

    if (len(zero_crossings) == 0):
        print("Uncontrolled Self-Energy result, no zero crossing. Returning infinity")
        print(loc_sigma)
        return float('inf')

    ind = zero_crossings[0]
    w_prev = w_points[ind]
    w_next = w_points[ind+1]
    sigma_prev = loc_sigma[ind]
    sigma_next = loc_sigma[ind+1]
    renorm = 1.0 - (sigma_next.imag - sigma_prev.imag) / (w_next - w_prev)
    print("w-cross: ", w_prev)

    return renorm

def contract(obj1, obj2):
    shape1 = obj1.data.shape
    shape2 = obj2.data.shape
    s1 = len(shape1)
    s2 = len(shape2)
    if shape1 == shape2:
        new_obj = obj1.copy()
        new_obj.data[:] = obj1.data * obj2.data
    elif s2 > s1:
        new_obj = obj1.copy()
        new_obj.data[:] = np.einsum('wkabcd,wkcd->wkab', obj2.data, obj1.data, optimize=True)
    elif s1 > s2:
        new_obj = obj2.copy()
        new_obj.data[:] = np.einsum('wkabcd,wkcd->wkab', obj1.data, obj2.data, optimize=True)
    else:
        raise ValueError("Contraction only implemented for Gf objects with same shape or differing by 2 indices")
    return new_obj # Returns a Gf object


def dot_t(diagram1, diagram2):
    shape1 = diagram1.obj_t.data.shape
    shape2 = diagram2.obj_t.data.shape
    s = len(shape1) - len(shape2)
    if s <= 0:
        new_obj = diagram1.copy()
    else:
        new_obj = diagram2.copy()
    new_obj.zero()
    if s == 0:
        new_obj = diagram1.copy()
        new_obj.obj_t.data[:] = diagram1.obj_t.data * diagram2.obj_t.data
        return new_obj
    elif abs(s) == 2:
        n_orbs = shape1[2] if s > 0 else shape2[2]
        for i in range(n_orbs):
            for j in range(n_orbs):
                if s > 0:
                    new_obj.obj_t.data[:] += diagram1.obj_t.data[:, i, j] * diagram2.obj_t.data[:]
                else:
                    new_obj.obj_t.data[:] += diagram1.obj_t.data[:] * diagram2.obj_t.data[:, i, j]
        return new_obj
    else:
        raise ValueError("Convolution Sum only implemented for diagrams differing by 0 or 2 indices")


def dot_tr(diagram1, diagram2):
    shape1 = diagram1.obj_tr.data.shape
    shape2 = diagram2.obj_tr.data.shape
    s = len(shape1) - len(shape2)
    if s <= 0:
        new_obj = diagram1.copy()
    else:
        new_obj = diagram2.copy()
    new_obj.zero()

    if s == 0:
        new_obj = diagram1.copy()
        new_obj.obj_tr.data[:] = diagram1.obj_tr.data * diagram2.obj_tr.data
        return new_obj
    elif abs(s) == 2:
        n_orbs = shape1[2] if s > 0 else shape2[2]
        for i in range(n_orbs):
            for j in range(n_orbs):
                if s > 0:
                    new_obj.obj_tr.data[:] += diagram1.obj_tr.data[:, :, i, j] * diagram2.obj_tr.data[:]
                else:
                    new_obj.obj_tr.data[:] += diagram1.obj_tr.data[:] * diagram2.obj_tr.data[:, :, i, j]
        return new_obj
    else:
        raise ValueError("Convolution Sum only implemented for diagrams differing by 0 or 2 indices")

def flip_k(diagram):
    if diagram.varspace not in ['wk', 'tr']:
        raise ValueError("flip_wk only works for diagrams in wk or tr space")

    # Create a copy
    G_flip = diagram.copy()

    # Get mesh dimensions
    mesh, _ = extract_mesh_and_bz(diagram.obj_wk)
    original_shape = diagram.obj_wk.data.shape

    # Reshape to mesh dimensions (nω, nkx, nky, nkz, orb1, orb2, ...)
    data = np.reshape(diagram.obj_wk.data, mesh + original_shape[2:])
    k_axes = tuple(range(1, len(mesh)))
    data = np.flip(data, axis=k_axes)

    # Reshape back to flat format
    G_flip.obj_wk.data[:] = np.reshape(data, original_shape)

    G_flip.wk_to_tr()

    return G_flip

def flip_w(diagram):
    if diagram.varspace not in ['wk', 'tr']:
        raise ValueError("flip_wk only works for diagrams in wk or tr space")

    # Create a copy
    G_flip = diagram.copy()

    # Get mesh dimensions
    mesh, _ = extract_mesh_and_bz(diagram.obj_wk)
    original_shape = diagram.obj_wk.data.shape

    # Reshape to mesh dimensions (nω, nkx, nky, nkz, orb1, orb2, ...)
    data = np.reshape(diagram.obj_wk.data, mesh + original_shape[2:])
    data = np.flip(data, axis=0)
    # Reshape back to flat format
    G_flip.obj_wk.data[:] = (G_flip.obj_wk.data[:] + np.reshape(data, original_shape)) / 2.0

    G_flip.wk_to_tr()

    return G_flip

def flip_wk(diagram):
    """
    Flip Green's function in momentum and frequency: G(k, iω) → G(-k, -iω).

    This operation is crucial for Eliashberg equations where Cooper pairs
    require G(k)G(-k) with opposite momenta.

    Args:
        diagram: Diagram object in wk space

    Returns:
        Diagram object with flipped data: G(-k, -iω)
    """
    if diagram.varspace not in ['wk', 'tr']:
        raise ValueError("flip_wk only works for diagrams in wk or tr space")

    # Create a copy
    G_flip = diagram.copy()

    # Get mesh dimensions
    mesh, _ = extract_mesh_and_bz(diagram.obj_wk)
    original_shape = diagram.obj_wk.data.shape

    # Reshape to mesh dimensions (nω, nkx, nky, nkz, orb1, orb2, ...)
    data = np.reshape(diagram.obj_wk.data, mesh + original_shape[2:])

    # Flip frequency: iω → -iω
    data = np.flip(data, axis=0)

    # Flip all k-axes: k → -k
    k_axes = tuple(range(1, len(mesh)))
    data = np.flip(data, axis=k_axes)

    # Reshape back to flat format
    G_flip.obj_wk.data[:] = np.reshape(data, original_shape)

    # Update tr representation
    G_flip.wk_to_tr()

    return G_flip


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
    elif isinstance(mesh, BrillouinZone):
        return "k"
    else:
        return "k"
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
        if isinstance(mesh, MeshBrillouinZone):
            raise ValueError("Mesh is only MeshBrillouinZone, missing frequency component")

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

    # Determine actual dimensionality from k_dims
    # For 2D systems, k_dims = (nkx, nky, 1), so trim BZ to 2x2
    # For 1D systems, k_dims = (nkx, 1, 1), so trim BZ to 1x1
    if k_dims[2] == 1 and k_dims[1] == 1:
        # 1D system
        brillouin_zone = brillouin_zone[:1, :1]
    elif k_dims[2] == 1:
        # 2D system
        brillouin_zone = brillouin_zone[:2, :2]
    # else: 3D system, keep full 3x3 matrix

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


def get_brillouin_zone():
    """
    Construct TRIQS BrillouinZone from config parameters.

    Returns:
        BrillouinZone: TRIQS BrillouinZone object based on cfg.brillouin_zone and cfg.dimension

    Example:
        >>> BZ = get_brillouin_zone()
        >>> k_mesh = MeshBrZone(BZ, n_k=60)
    """
    # Get reciprocal lattice vectors from config
    BZ_vectors = np.array(cfg.brillouin_zone)

    # Extract only the relevant dimensions based on cfg.dimension
    if cfg.dimension == 2:
        # For 2D systems, use only the first 2x2 block
        BZ_vectors = BZ_vectors[:2, :2]

    # Compute real-space lattice vectors from reciprocal vectors
    # Using the relation: a_i · b_j = 2π δ_ij
    a = 2 * np.pi * np.linalg.inv(BZ_vectors.T).T

    # Create BravaisLattice with appropriate dimensionality
    if cfg.dimension == 2:
        # 2D system
        bl = BravaisLattice(units=[[a[0,0], a[0,1]],
                                   [a[1,0], a[1,1]]])
    else:
        # 3D system
        bl = BravaisLattice(units=[[a[0,0], a[0,1], a[0,2]],
                                   [a[1,0], a[1,1], a[1,2]],
                                   [a[2,0], a[2,1], a[2,2]]])

    return BrillouinZone(bl)


