import numpy as np
import spglib
from scipy.spatial import ConvexHull

def get_ibz(lattice, positions, numbers, mesh=(10, 10, 10)):
    # Define crystal structure
    cell = (lattice, positions, numbers)
    
    # Compute reciprocal lattice vectors
    volume = np.dot(lattice[0], np.cross(lattice[1], lattice[2]))
    reciprocal_lattice = 2 * np.pi * np.linalg.inv(lattice).T  # Convert to reciprocal basis
    
    # Get symmetry operations using spglib
    symmetry = spglib.get_symmetry(cell)
    rotations = np.array(symmetry['rotations'])
    translations = np.array(symmetry['translations'])
    
    # Generate k-points in the full Brillouin zone
    _, mesh_points = spglib.get_ir_reciprocal_mesh(mesh, cell)  # Extract correct k-points
    mesh_points = np.array(mesh_points, dtype=float)
    
    # Convert k-points to absolute reciprocal lattice units
    mesh_points = np.dot(mesh_points, reciprocal_lattice)
    
    # Wrap k-points into the first Brillouin zone (limit within -π/a to π/a)
    bz_limit = np.linalg.norm(reciprocal_lattice, axis=1) / 2
    mesh_points = np.mod(mesh_points + bz_limit, 2 * bz_limit) - bz_limit
    
    # Map each k-point to its irreducible representation
    unique_kpoints = set()
    for k in mesh_points:
        k_transformed = [tuple(np.round(k, decimals=8))]  # Ensure tuple conversion and avoid floating point errors
        for R, t in zip(rotations, translations):
            k_sym = np.dot(R, k) + t  # Apply symmetry operation
            k_sym = np.mod(k_sym + bz_limit, 2 * bz_limit) - bz_limit  # Wrap into first BZ
            k_transformed.append(tuple(np.round(k_sym, decimals=8)))
        
        # Keep only one representative per symmetry-equivalent group
        unique_kpoints.add(min(k_transformed))
    
    unique_kpoints = np.array(list(unique_kpoints))
    
    # Select three independent vectors from the IBZ to form a 3x3 basis
    if len(unique_kpoints) >= 3:
        ibz_basis_vectors = unique_kpoints[:3]  # Select the first three k-points as basis vectors
    else:
        ibz_basis_vectors = unique_kpoints  # Use what is available if fewer than 3 exist
    
    return unique_kpoints, ibz_basis_vectors, reciprocal_lattice

# Example usage (Simple cubic lattice)
lattice = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
positions = np.array([[0, 0, 0]])
numbers = [1]

ibz_kpoints, ibz_basis_vectors, reciprocal_lattice = get_ibz(lattice, positions, numbers)
print("Irreducible k-points:", ibz_kpoints)
print("IBZ Basis Vectors (3x3):", ibz_basis_vectors)
print("Reciprocal Lattice Vectors:", reciprocal_lattice)
