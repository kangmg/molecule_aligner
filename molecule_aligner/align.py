"""
Core alignment and merging functionality for molecular trajectories.
"""

from typing import List, Dict, Union, Optional
from pathlib import Path
import ase
from ase import Atoms
from ase.io import read, write

from .utils import align_atoms_to_reference


def load_trajectory(traj_dict: Dict[str, Union[str, List[Atoms], bool]]) -> List[Atoms]:
    """
    Load trajectory from dictionary specification.
    
    Args:
        traj_dict: Dictionary containing either 'traj_path' or 'traj' key
        
    Returns:
        List of Atoms objects representing the trajectory
        
    Raises:
        ValueError: If neither 'traj_path' nor 'traj' is provided, or both are provided
    """
    has_path = 'traj_path' in traj_dict and traj_dict['traj_path'] is not None
    has_traj = 'traj' in traj_dict and traj_dict['traj'] is not None
    
    if has_path and has_traj:
        raise ValueError("Both 'traj_path' and 'traj' provided. Only one should be specified.")
    
    if not has_path and not has_traj:
        raise ValueError("Neither 'traj_path' nor 'traj' provided. One must be specified.")
    
    if has_path:
        # Load trajectory from file
        traj_path = traj_dict['traj_path']
        trajectory = read(traj_path, ':')
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
    else:
        # Use provided trajectory
        trajectory = traj_dict['traj']
    
    # Reverse trajectory if requested
    if traj_dict.get('reverse', False):
        trajectory = trajectory[::-1]
    
    return trajectory


def align_trajectory_to_reference(
    trajectory: List[Atoms], 
    reference_frame: Atoms,
    base_indices: List[int]
) -> List[Atoms]:
    """
    Align all frames in a trajectory to a reference frame using base indices.
    
    Args:
        trajectory: List of Atoms objects to align
        reference_frame: Reference Atoms object for alignment
        base_indices: List of atom indices to use for alignment
        
    Returns:
        List of aligned Atoms objects
    """
    aligned_trajectory = []
    
    for frame in trajectory:
        # Validate that base_indices exist in both frames
        if max(base_indices) >= len(frame):
            raise ValueError(f"Base indices {base_indices} exceed frame size {len(frame)}")
        if max(base_indices) >= len(reference_frame):
            raise ValueError(f"Base indices {base_indices} exceed reference frame size {len(reference_frame)}")
        
        aligned_frame = align_atoms_to_reference(
            atoms=frame,
            reference=reference_frame,
            atoms_indices=base_indices,
            reference_indices=base_indices
        )
        aligned_trajectory.append(aligned_frame)
    
    return aligned_trajectory


def align_and_merge_reactions(
    reactions: List[Dict[str, Union[str, List[Atoms], bool, List[int]]]],
    output_path: Optional[str] = None,
    reference: str = 'first'
) -> List[Atoms]:
    """
    Align and merge multiple reaction trajectories based on shared base molecule.
    
    Args:
        reactions: List of dictionaries, each containing:
            - 'traj_path' (optional): Path to .traj file
            - 'traj' (optional): List of ASE Atoms objects  
            - 'reverse' (bool): Whether to reverse the trajectory
            - 'base_indices' (List[int]): Atom indices for alignment
        output_path: Optional path to save merged trajectory as .extxyz
        reference: Reference frame selection ('first' or 'reactant')
        
    Returns:
        Merged trajectory as list of Atoms objects
        
    Raises:
        ValueError: If reactions list is empty or contains invalid data
    """
    if not reactions:
        raise ValueError("No reactions provided")
    
    # Load all trajectories
    trajectories = []
    base_indices_list = []
    
    for i, reaction_dict in enumerate(reactions):
        if 'base_indices' not in reaction_dict:
            raise ValueError(f"Reaction {i} missing required 'base_indices' field")
        
        trajectory = load_trajectory(reaction_dict)
        trajectories.append(trajectory)
        base_indices_list.append(reaction_dict['base_indices'])
    
    if not trajectories:
        raise ValueError("No valid trajectories loaded")
    
    # Determine reference frame
    if reference == 'first':
        reference_frame = trajectories[0][0]
        reference_base_indices = base_indices_list[0]
    elif reference == 'reactant':
        # Use first frame of first trajectory as reactant
        reference_frame = trajectories[0][0]  
        reference_base_indices = base_indices_list[0]
    else:
        raise ValueError(f"Invalid reference option '{reference}'. Use 'first' or 'reactant'.")
    
    # Align all trajectories to the reference
    merged_trajectory = []
    
    for i, (trajectory, base_indices) in enumerate(zip(trajectories, base_indices_list)):
        # Validate base indices consistency
        if len(base_indices) != len(reference_base_indices):
            raise ValueError(
                f"Trajectory {i} base_indices length ({len(base_indices)}) "
                f"does not match reference ({len(reference_base_indices)})"
            )
        
        aligned_trajectory = align_trajectory_to_reference(
            trajectory=trajectory,
            reference_frame=reference_frame,
            base_indices=base_indices
        )
        
        merged_trajectory.extend(aligned_trajectory)
    
    # Save to file if output path is provided
    if output_path:
        output_path = Path(output_path)
        # Ensure .extxyz extension
        if output_path.suffix.lower() != '.extxyz':
            output_path = output_path.with_suffix('.extxyz')
        
        write(str(output_path), merged_trajectory)
    
    return merged_trajectory