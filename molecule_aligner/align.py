"""
Core alignment and merging functionality for molecular trajectories.
Enhanced with single frame interpolation support.
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


def detect_input_type(reaction_dict: Dict) -> str:
    """
    Detect input type: 'trajectory', 'interpolation'.
    
    Args:
        reaction_dict: Input dictionary to analyze
        
    Returns:
        'trajectory': Multi-frame trajectory (existing behavior)
        'interpolation': Two structures for interpolation
        
    Raises:
        ValueError: If input format is invalid or ambiguous
    """
    
    # Check for trajectory inputs
    has_traj_path = 'traj_path' in reaction_dict and reaction_dict['traj_path'] is not None
    has_traj = 'traj' in reaction_dict and reaction_dict['traj'] is not None
    has_trajectory = has_traj_path or has_traj
    
    # Check for interpolation inputs
    has_start_frame = 'start_frame' in reaction_dict or 'start_frame_path' in reaction_dict
    has_end_frame = 'end_frame' in reaction_dict or 'end_frame_path' in reaction_dict
    has_interpolation = has_start_frame and has_end_frame
    
    if has_trajectory and has_interpolation:
        raise ValueError(
            "Cannot specify both trajectory and interpolation inputs. "
            "Use either ('traj'/'traj_path') OR ('start_frame'+'end_frame')."
        )
    
    if has_interpolation:
        return 'interpolation'
    elif has_trajectory:
        return 'trajectory'
    else:
        raise ValueError(
            "Must specify either trajectory inputs ('traj'/'traj_path') "
            "or interpolation inputs ('start_frame'+'end_frame')."
        )


def load_trajectory_enhanced(reaction_dict: Dict) -> List[Atoms]:
    """
    Enhanced trajectory loading supporting both trajectories and interpolation.
    
    Args:
        reaction_dict: Dictionary with trajectory or interpolation specification
        
    Returns:
        List of Atoms objects representing the trajectory
        
    Raises:
        ValueError: If input format is invalid
    """
    
    input_type = detect_input_type(reaction_dict)
    
    if input_type == 'trajectory':
        # Use existing load_trajectory function for backward compatibility
        return load_trajectory(reaction_dict)
    
    elif input_type == 'interpolation':
        from .interpolate import create_interpolated_trajectory
        
        # Load start frame
        if 'start_frame' in reaction_dict:
            start_atoms = reaction_dict['start_frame']
            if not isinstance(start_atoms, Atoms):
                raise ValueError("start_frame must be an ASE Atoms object")
        elif 'start_frame_path' in reaction_dict:
            start_atoms = read(reaction_dict['start_frame_path'])
        else:
            raise ValueError("Must specify either 'start_frame' or 'start_frame_path'")
        
        # Load end frame  
        if 'end_frame' in reaction_dict:
            end_atoms = reaction_dict['end_frame']
            if not isinstance(end_atoms, Atoms):
                raise ValueError("end_frame must be an ASE Atoms object")
        elif 'end_frame_path' in reaction_dict:
            end_atoms = read(reaction_dict['end_frame_path'])
        else:
            raise ValueError("Must specify either 'end_frame' or 'end_frame_path'")
        
        # Get interpolation parameters
        n_frames = reaction_dict.get('n_frames', 10)
        method = reaction_dict.get('interpolation_method', 'idpp')
        base_indices = reaction_dict.get('base_indices')
        options = reaction_dict.get('interpolation_options', {})
        
        # Validate parameters
        if n_frames < 2:
            raise ValueError("n_frames must be at least 2")
        
        # Create interpolated trajectory
        trajectory = create_interpolated_trajectory(
            start_atoms=start_atoms,
            end_atoms=end_atoms,
            n_frames=n_frames,
            method=method,
            base_indices=base_indices,
            **options
        )
        
        # Apply reverse if requested
        if reaction_dict.get('reverse', False):
            trajectory = trajectory[::-1]
        
        return trajectory
    
    else:
        raise ValueError(f"Unknown input type: {input_type}")


def align_and_merge_reactions_enhanced(
    reactions: List[Dict[str, Union[str, List[Atoms], bool, List[int]]]],
    output_path: Optional[str] = None,
    reference: str = 'first'
) -> List[Atoms]:
    """
    Enhanced version of align_and_merge_reactions supporting both trajectories and interpolation.
    
    This maintains full backward compatibility while adding interpolation support.
    
    Args:
        reactions: List of dictionaries, each containing either:
            TRAJECTORY FORMAT:
            - 'traj_path' (optional): Path to .traj file
            - 'traj' (optional): List of ASE Atoms objects  
            - 'reverse' (bool): Whether to reverse the trajectory
            - 'base_indices' (List[int]): Atom indices for alignment
            
            INTERPOLATION FORMAT:
            - 'start_frame' OR 'start_frame_path': Starting structure
            - 'end_frame' OR 'end_frame_path': Ending structure
            - 'n_frames' (int): Number of interpolated frames (default: 10)
            - 'interpolation_method' (str): 'linear' or 'idpp' (default: 'idpp')
            - 'interpolation_options' (dict): Method-specific options
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
    
    # Load all trajectories (now supports interpolation)
    trajectories = []
    base_indices_list = []
    
    for i, reaction_dict in enumerate(reactions):
        if 'base_indices' not in reaction_dict:
            raise ValueError(f"Reaction {i} missing required 'base_indices' field")
        
        # Use enhanced loader that supports both trajectories and interpolation
        trajectory = load_trajectory_enhanced(reaction_dict)
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


# Keep original function for backward compatibility
def align_and_merge_reactions(
    reactions: List[Dict[str, Union[str, List[Atoms], bool, List[int]]]],
    output_path: Optional[str] = None,
    reference: str = 'first'
) -> List[Atoms]:
    """
    Original align_and_merge_reactions function - maintained for backward compatibility.
    
    For new features including interpolation, use align_and_merge_reactions_enhanced().
    
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
    """
    
    if not reactions:
        raise ValueError("No reactions provided")
    
    # Load all trajectories using original method
    trajectories = []
    base_indices_list = []
    
    for i, reaction_dict in enumerate(reactions):
        if 'base_indices' not in reaction_dict:
            raise ValueError(f"Reaction {i} missing required 'base_indices' field")
        
        trajectory = load_trajectory(reaction_dict)  # Original loader only
        trajectories.append(trajectory)
        base_indices_list.append(reaction_dict['base_indices'])
    
    if not trajectories:
        raise ValueError("No valid trajectories loaded")
    
    # Determine reference frame
    if reference == 'first':
        reference_frame = trajectories[0][0]
        reference_base_indices = base_indices_list[0]
    elif reference == 'reactant':
        reference_frame = trajectories[0][0]  
        reference_base_indices = base_indices_list[0]
    else:
        raise ValueError(f"Invalid reference option '{reference}'. Use 'first' or 'reactant'.")
    
    # Align all trajectories to the reference
    merged_trajectory = []
    
    for i, (trajectory, base_indices) in enumerate(zip(trajectories, base_indices_list)):
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
        if output_path.suffix.lower() != '.extxyz':
            output_path = output_path.with_suffix('.extxyz')
        write(str(output_path), merged_trajectory)
    
    return merged_trajectory


def create_reaction_path(
    start_structure: Union[Atoms, str],
    end_structure: Union[Atoms, str],
    base_indices: List[int],
    n_frames: int = 20,
    method: str = 'idpp',
    output_path: Optional[str] = None,
    **method_options
) -> List[Atoms]:
    """
    Simplified interface for creating reaction paths between two structures.
    
    Args:
        start_structure: Starting structure (Atoms object or file path)
        end_structure: Ending structure (Atoms object or file path)  
        base_indices: Atom indices for alignment
        n_frames: Number of frames to generate (default: 20)
        method: Interpolation method 'linear' or 'idpp' (default: 'idpp')
        output_path: Optional output file path
        **method_options: Method-specific options:
            For IDPP: fmax, steps, optimizer, mic
            For linear: mic, interpolate_cell, use_scaled_coord
        
    Returns:
        Interpolated reaction path as list of Atoms objects
        
    Raises:
        ValueError: If input parameters are invalid
    """
    
    from .interpolate import create_interpolated_trajectory
    
    # Load structures
    if isinstance(start_structure, str):
        start_atoms = read(start_structure)
    elif isinstance(start_structure, Atoms):
        start_atoms = start_structure
    else:
        raise ValueError("start_structure must be Atoms object or file path string")
        
    if isinstance(end_structure, str):
        end_atoms = read(end_structure)
    elif isinstance(end_structure, Atoms):
        end_atoms = end_structure
    else:
        raise ValueError("end_structure must be Atoms object or file path string")
    
    # Validate base_indices
    if not isinstance(base_indices, list) or len(base_indices) == 0:
        raise ValueError("base_indices must be a non-empty list")
    
    # Create interpolated trajectory
    trajectory = create_interpolated_trajectory(
        start_atoms=start_atoms,
        end_atoms=end_atoms,
        n_frames=n_frames,
        method=method,
        base_indices=base_indices,
        **method_options
    )
    
    # Save if requested
    if output_path:
        write(output_path, trajectory)
    
    return trajectory