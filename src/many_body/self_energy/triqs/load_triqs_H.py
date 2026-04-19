from triqs.gf import *
from triqs.gf.meshes import MeshDLRImFreq
from triqs_tprf.tight_binding import TBLattice

import firefly.config as cfg

hamiltonian_type = cfg.hamiltonian

def create_dlr_meshes(e_k, beta, statistic='Fermion'):
    """
    Create DLR frequency mesh from energy dispersion.

    Args:
        e_k: Energy dispersion from get_energy_mesh()
        beta: Inverse temperature
        statistic: 'Fermion' or 'Boson'

    Returns:
        MeshDLRImFreq with wmax = 1.2 * (max_ek - min_ek), eps=1e-14
    """
    emax = e_k.data.max().real
    emin = e_k.data.min().real
    w_max = 1.2 * (emax - emin)
    eps = 1e-14

    mesh = MeshDLRImFreq(beta=beta, statistic=statistic, w_max=w_max, eps=eps)
    ws = [w.value for w in mesh]
    print("Max w: ", max(ws))
    print("Min w: ", min(ws))
    return mesh


def get_energy_mesh():
    if hamiltonian_type == "tight_binding":
        return get_TB()
    if hamiltonian_type == "emery":
        return get_emery()
    else:
        raise ValueError(f"Hamiltonian type {hamiltonian_type} not recognized.")

def get_emery():
    """
    Create the Emery model (3-band model for cuprate superconductors).

    The model has 3 orbitals per unit cell:
    - Orbital 0: Cu d_{x²-y²} at (0,0)
    - Orbital 1: O p_x at (0.5, 0)
    - Orbital 2: O p_y at (0, 0.5)

    Standard parameters (in eV):
    - ε_d = 0.0 : Cu d-orbital on-site energy (reference)
    - ε_p = 3.6 : O p-orbital on-site energy
    - t_pd = 1.3 : Cu-O hopping
    - t_pp = 0.65 : O-O hopping

    Returns:
        tuple: (H_r, kmesh, e_k)
    """
    import firefly.config as cfg

    # Standard Emery model parameters
    eps_d = 0.0    # Cu d-orbital energy (reference)
    eps_p = 3.6    # O p-orbital energy
    t_pd = 1.3     # Cu d - O p hopping
    t_pp = 0.65    # O p - O p hopping

    Nk = cfg.k_mesh[0]

    # Create 3-orbital tight binding Hamiltonian on 2D square lattice
    # Unit cell contains: [d, px, py] orbitals
    H_r = TBLattice(
        units=[
            (1, 0, 0),  # a₁ lattice vector
            (0, 1, 0),  # a₂ lattice vector
        ],
        orbital_positions=[
            (0, 0, 0),      # Cu d at origin
            (0.5, 0, 0),    # O px along x
            (0, 0.5, 0),    # O py along y
        ],
        orbital_names=['d', 'px', 'py'],
        hoppings={
            # On-site energies
            (0, 0): [[eps_d, 0, 0],
                     [0, eps_p, 0],
                     [0, 0, eps_p]],

            # Cu d - O px hopping (±x direction)
            # d to px in same unit cell
            # Factor includes phase from orbital positions
            (+1, 0): [[0, -t_pd, 0],
                      [-t_pd, 0, 0],
                      [0, 0, 0]],
            (-1, 0): [[0, -t_pd, 0],
                      [-t_pd, 0, 0],
                      [0, 0, 0]],

            # Cu d - O py hopping (±y direction)
            (0, +1): [[0, 0, -t_pd],
                      [0, 0, 0],
                      [-t_pd, 0, 0]],
            (0, -1): [[0, 0, -t_pd],
                      [0, 0, 0],
                      [-t_pd, 0, 0]],

            # O px - O py hopping (diagonal)
            (+1, +1): [[0, 0, 0],
                       [0, 0, t_pp],
                       [0, t_pp, 0]],
            (+1, -1): [[0, 0, 0],
                       [0, 0, -t_pp],
                       [0, -t_pp, 0]],
            (-1, +1): [[0, 0, 0],
                       [0, 0, -t_pp],
                       [0, -t_pp, 0]],
            (-1, -1): [[0, 0, 0],
                       [0, 0, t_pp],
                       [0, t_pp, 0]],
        })

    # Create k-mesh and compute dispersion
    kmesh = H_r.get_kmesh(n_k=Nk)
    e_k = H_r.fourier(kmesh)

    return H_r, kmesh, e_k

def get_TB():
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

