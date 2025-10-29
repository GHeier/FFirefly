from triqs.gf import *
from triqs_tprf.tight_binding import TBLattice


def get_energy_mesh():
    """
    Create tight-binding lattice and energy mesh based on firefly.config parameters.

    Returns:
        tuple: (H_r, kmesh, e_k) where
            - H_r: TBLattice object (2D or 3D based on cfg.dimension)
            - kmesh: k-space mesh
            - e_k: energy dispersion on k-mesh
    """
    # Import config inside function to avoid circular import
    import firefly.config as cfg

    # Parameters
    t = cfg.t0
    t1 = cfg.t1
    t2 = cfg.t2
    Nk = cfg.k_mesh[0]
    dim = cfg.dimension

    if dim == 2:
        # Create 2D tight binding Hamiltonian
        H_r = TBLattice(
            units=[
                (1, 0, 0),  # basis vector in the x-direction
                (0, 1, 0),  # basis vector in the y-direction
            ],
            hoppings={
                (+1, 0): [[-t]],   # nearest-neighbor hopping in +x
                (-1, 0): [[-t]],   # nearest-neighbor hopping in -x
                (0, +1): [[-t]],   # nearest-neighbor hopping in +y
                (0, -1): [[-t]],   # nearest-neighbor hopping in -y
                (+1, +1): [[-t1]], # next-nearest-neighbor hopping
                (+1, -1): [[-t1]],
                (-1, +1): [[-t1]],
                (-1, -1): [[-t1]],
                (+2, 0): [[-t2]],  # next-next-nearest-neighbor hopping
                (-2, 0): [[-t2]],
                (0, +2): [[-t2]],
                (0, -2): [[-t2]],
            })

    elif dim == 3:
        # Create 3D tight binding Hamiltonian
        H_r = TBLattice(
            units=[
                (1, 0, 0),
                (0, 1, 0),
                (0, 0, 1),
            ],
            hoppings={
                (+1, 0, 0): [[-t]],  # nearest-neighbor hopping
                (-1, 0, 0): [[-t]],
                (0, +1, 0): [[-t]],
                (0, -1, 0): [[-t]],
                (0, 0, +1): [[-t]],
                (0, 0, -1): [[-t]],
            })

    else:
        raise ValueError(f"Unsupported dimension: {dim}. Only 2D and 3D are supported.")

    # Create k-mesh and compute dispersion
    kmesh = H_r.get_kmesh(n_k=Nk)
    e_k = H_r.fourier(kmesh)

    return H_r, kmesh, e_k

