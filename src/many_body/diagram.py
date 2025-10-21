from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
from triqs.gf import MeshProduct, MeshBrillouinZone
from triqs_tprf.lattice import fourier_tr_to_wr, fourier_wk_to_wr, fourier_wr_to_tr, fourier_wr_to_wk, chi_wr_from_chi_tr, chi_wk_from_chi_wr, chi_tr_from_chi_wr, chi_wr_from_chi_wk
import numpy as np

import firefly as fly
import firefly.config as cfg

class Diagram:
    def __init__(self, obj, statistic):
        varspace = describe_mesh(obj)
        if varspace not in ['wk', 'tr']:
            print(f"Diagram initialized in {varspace} space.")
            raise ValueError("Object wk-mesh must be MeshDLRImFreq or MeshDLRImTime")

        self.statistic = statistic
        self.shape = obj.data.shape
        self.nw = self.shape[0]
        self.nk = self.shape[1]
        self.ind_dim = len(self.shape) - 2

        # Extract w-points from mesh (Matsubara frequencies)
        mesh_w = obj.mesh.components[0]
        # For Matsubara frequencies, use imaginary part
        self.w_points = np.array([float(iw.imag) for iw in mesh_w], dtype=np.float32)

        if varspace == 'wk':
            self.obj_wk = obj
            self.wk_to_tr()
        elif varspace == 'tr':
            self.obj_tr = obj
            self.tr_to_wk()
        else:
            raise ValueError("varspace must be 'wk' or 'tr'")

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
    def save_w(self, filename):
        print(self.obj_wk.data.shape)
        obj_w = np.sum(self.obj_wk.data, axis=1) / self.nk  # Sum over k-points
        print(obj_w.shape)
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

    #elif isinstance(mesh, MeshDLRImFreq):
    #    return "w"
    #elif isinstance(mesh, MeshDLRImTime):
    #    return "t"
    #elif isinstance(mesh, BrillouinZoneMesh):
    #    return "k"
    else:
        return "unknown mesh type"

