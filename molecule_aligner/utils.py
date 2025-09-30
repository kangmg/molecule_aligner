"""
Utility functions for molecular alignment.
"""

from typing import List
import numpy as np
from ase import Atoms


def kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """
    Kabsch algorithm to compute the optimal rotation matrix that aligns two sets of points P and Q.
    
    Args:
        P: First set of points (shape: N x 3)
        Q: Second set of points (shape: N x 3)
        
    Returns:
        Optimal rotation matrix (3 x 3)
    """
    # Center the points
    P_centered = P - np.mean(P, axis=0)
    Q_centered = Q - np.mean(Q, axis=0)

    # Compute the covariance matrix
    C = np.dot(np.transpose(P_centered), Q_centered)
    
    # Singular value decomposition
    V, S, W = np.linalg.svd(C)
    
    # Correct for reflection
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    
    # Compute rotation matrix
    U = np.dot(V, W)
    return U


def align_atoms_to_reference(
    atoms: Atoms, 
    reference: Atoms, 
    atoms_indices: List[int], 
    reference_indices: List[int]
) -> Atoms:
    """
    Align atoms object to a reference structure using specified atom indices.
    
    Args:
        atoms: Atoms object to be aligned
        reference: Reference atoms object
        atoms_indices: List of atom indices in atoms to use for alignment
        reference_indices: List of atom indices in reference to use for alignment
        
    Returns:
        Aligned copy of the atoms object
    """
    aligned_atoms = atoms.copy()
    
    # Get coordinates of alignment atoms
    atoms_coords = atoms.get_positions()[atoms_indices]
    ref_coords = reference.get_positions()[reference_indices]
    
    # Calculate centroids
    atoms_centroid = np.mean(atoms_coords, axis=0)
    ref_centroid = np.mean(ref_coords, axis=0)
    
    # Center the atoms
    aligned_atoms.positions -= atoms_centroid
    
    # Calculate rotation matrix using centered coordinates
    atoms_coords_centered = atoms_coords - atoms_centroid
    ref_coords_centered = ref_coords - ref_centroid
    
    rotation_matrix = kabsch(atoms_coords_centered, ref_coords_centered)
    
    # Apply rotation
    aligned_atoms.positions = np.dot(aligned_atoms.positions, rotation_matrix)
    
    # Translate to reference position
    aligned_atoms.positions += ref_centroid
    
    return aligned_atoms