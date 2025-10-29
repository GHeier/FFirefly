import numpy as np
from triqs.gf import Gf
from triqs.gf.meshes import MeshDLRImFreq
from triqs.gf import MeshProduct, MeshBrillouinZone
from triqs.lattice.lattice_tools import BrillouinZone, BravaisLattice

def fill_triqs_from_field(Gf_obj, field):
    print("Filling TRIQS Green's function from FFirefly field...")
    """
    Fill a TRIQS Green's function with interpolated values from a FFirefly field.

    Args:
        Gf_obj: TRIQS Green's function object to fill
        field: FFirefly Field object to read from
    """
    mesh = Gf_obj.mesh
    # Check if mesh is DLR imfreq only, or a product mesh with DLR imfreq and BZ mesh
    if isinstance(mesh, MeshProduct):
        # Product mesh: extract frequency and k-mesh components
        mesh_iw, mesh_k = mesh.components
        if not isinstance(mesh_iw, MeshDLRImFreq):
            raise ValueError("Expected first mesh component to be MeshDLRImFreq")
        if not isinstance(mesh_k, MeshBrillouinZone):
            raise ValueError("Expected second mesh component to be MeshBrillouinZone")
        BZ_from_mesh = mesh_k.bz.units
        BZ = np.array(BZ_from_mesh)
    elif isinstance(mesh, MeshDLRImFreq):
        # Simple frequency mesh only
        mesh_iw = mesh
        mesh_k = None
    else:
        raise ValueError(f"Unsupported mesh type: {type(mesh)}")

    if mesh_k is None:
        for iw_idx, iw in enumerate(mesh_iw):
            w_val = float(iw.real)
            values = field(w_val)
            if Gf_obj.target_shape[0] > 1:
                Gf_obj.data[iw_idx, :, :] = values
            else:
                Gf_obj.data[iw_idx] = values
    else:

        print("Using k-mesh for filling TRIQS Gf...")
        # Get all k-points from TRIQS mesh as a single Nx3 numpy array
        print("mesh_k ", mesh_k)
        nx, ny, nz = mesh_k.dims
        kx, ky, kz = np.meshgrid(np.arange(nx) / nx,
                                np.arange(ny) / ny,
                                np.arange(nz) / nz,
                                indexing='ij')
        kpts = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], axis=-1)
        k_points = kpts @ BZ.T
        nk = len(k_points)

        print(f"Total k-points: nk={nk}")
        # Check if field is matrix-valued
        is_matrix = hasattr(field, 'mat_dim') and field.mat_dim > 1
        print(f"Filling TRIQS Gf: nk={nk}, is_matrix={is_matrix}")

        # Loop over frequencies and evaluate all k-points at once
        for iw_idx, iw in enumerate(mesh_iw):
            print(iw_idx)
            # Get frequency value (imaginary part for Matsubara)
            w_val = float(iw.real)

            # Batch evaluate field at all k-points for this frequency
            # field(k_points, w) returns array of values for all k-points
            print("finding vals")
            values = field(k_points, w_val)
            print("filling obj")

            # Assign to Green's function data
            if is_matrix:
                # Matrix-valued field: values shape is (nk, mat_dim, mat_dim)
                Gf_obj.data[iw_idx, :, :, :] = values
            else:
                # Scalar field: values shape is (nk,)
                # Need to handle different target_shape dimensions
                target_ndim = len(Gf_obj.target_shape)
                if target_ndim == 0:
                    # Pure scalar, shape (nk,)
                    Gf_obj.data[iw_idx, :] = values
                elif target_ndim == 2:
                    # 2-index tensor, shape (nk, 1, 1)
                    Gf_obj.data[iw_idx, :, 0, 0] = values
                elif target_ndim == 4:
                    # 4-index tensor (chi), shape (nk, 1, 1, 1, 1)
                    Gf_obj.data[iw_idx, :, 0, 0, 0, 0] = values
                else:
                    raise ValueError(f"Unsupported target_shape dimension: {target_ndim}")



