"""
Interpolation utilities leveraging ASE methods for single frame interpolation.
"""

from typing import List, Optional, Union, Dict, Any
import numpy as np
from ase import Atoms
from ase.mep.neb import interpolate, idpp_interpolate
from ase.optimize import MDMin, LBFGS, FIRE
from .utils import align_atoms_to_reference


def create_interpolated_trajectory(
    start_atoms: Atoms,
    end_atoms: Atoms,
    n_frames: int,
    method: str = 'idpp',
    base_indices: Optional[List[int]] = None,
    **method_options
) -> List[Atoms]:
    """
    Create interpolated trajectory between two structures using ASE methods.
    
    Args:
        start_atoms: Initial structure
        end_atoms: Final structure  
        n_frames: Total frames (including start/end)
        method: 'linear' or 'idpp'
        base_indices: Atoms to align before interpolation
        **method_options: Method-specific parameters
            For IDPP: fmax, steps, optimizer, mic
            For linear: mic, interpolate_cell, use_scaled_coord
        
    Returns:
        List of interpolated Atoms objects
        
    Raises:
        ValueError: If method is not supported or n_frames < 2
    """
    
    if n_frames < 2:
        raise ValueError("n_frames must be at least 2 (start and end)")
    
    # Step 1: Align structures if base_indices provided
    if base_indices is not None:
        if len(base_indices) == 0:
            raise ValueError("base_indices cannot be empty")
        
        # Validate indices
        max_start_idx = len(start_atoms) - 1
        max_end_idx = len(end_atoms) - 1
        
        if max(base_indices) > max_start_idx:
            raise ValueError(f"base_indices {max(base_indices)} exceeds start_atoms size {len(start_atoms)}")
        if max(base_indices) > max_end_idx:
            raise ValueError(f"base_indices {max(base_indices)} exceeds end_atoms size {len(end_atoms)}")
        
        end_aligned = align_atoms_to_reference(
            atoms=end_atoms,
            reference=start_atoms,
            atoms_indices=base_indices,
            reference_indices=base_indices
        )
    else:
        end_aligned = end_atoms.copy()
    
    # Step 2: Create image list for interpolation
    images = [start_atoms.copy()]
    
    # Add intermediate blank images
    for i in range(1, n_frames - 1):
        intermediate = start_atoms.copy()
        images.append(intermediate)
    
    images.append(end_aligned.copy())
    
    # Step 3: Apply interpolation method
    if method == 'linear':
        # Use ASE linear interpolation
        mic = method_options.get('mic', False)
        interpolate_cell = method_options.get('interpolate_cell', False)
        use_scaled_coord = method_options.get('use_scaled_coord', False)
        apply_constraint = method_options.get('apply_constraint', None)
        
        interpolate(
            images, 
            mic=mic,
            interpolate_cell=interpolate_cell,
            use_scaled_coord=use_scaled_coord,
            apply_constraint=apply_constraint
        )
        
    elif method == 'idpp':
        # Use ASE IDPP interpolation
        fmax = method_options.get('fmax', 0.1)
        steps = method_options.get('steps', 100)
        mic = method_options.get('mic', False)
        optimizer_name = method_options.get('optimizer', 'MDMin')
        
        # Get optimizer class
        optimizer_map = {
            'MDMin': MDMin,
            'LBFGS': LBFGS, 
            'FIRE': FIRE
        }
        
        if optimizer_name in optimizer_map:
            optimizer = optimizer_map[optimizer_name]
        else:
            available = list(optimizer_map.keys())
            raise ValueError(f"Unsupported optimizer: {optimizer_name}. Available: {available}")
        
        try:
            # Apply IDPP interpolation
            idpp_interpolate(
                images,
                fmax=fmax,
                optimizer=optimizer,
                mic=mic,
                steps=steps
            )
        except Exception as e:
            # Fallback to linear interpolation if IDPP fails
            print(f"Warning: IDPP interpolation failed ({e}), falling back to linear interpolation")
            interpolate(images, mic=mic)
        
    else:
        available_methods = ['linear', 'idpp']
        raise ValueError(f"Unsupported method: {method}. Available: {available_methods}")
    
    return images


def interpolation_quality_check(
    trajectory: List[Atoms],
    start_atoms: Optional[Atoms] = None,
    end_atoms: Optional[Atoms] = None,
    method_name: str = "Unknown"
) -> Dict[str, float]:
    """
    Assess quality of interpolated trajectory.
    
    Args:
        trajectory: Interpolated trajectory to assess
        start_atoms: Original start structure (for accuracy check)
        end_atoms: Original end structure (for accuracy check)
        method_name: Name of interpolation method used
        
    Returns:
        Dictionary with quality metrics:
        - n_frames: Number of frames in trajectory
        - start_rmsd: RMSD between trajectory start and original start
        - end_rmsd: RMSD between trajectory end and original end  
        - max_displacement: Maximum single-step displacement
        - total_path_length: Total path length in coordinate space
        - smoothness_score: Path smoothness metric (lower = smoother)
        - avg_step_size: Average step size between frames
    """
    
    if not trajectory:
        return {'error': 'Empty trajectory'}
    
    metrics = {
        'method': method_name,
        'n_frames': len(trajectory),
        'start_rmsd': 0.0,
        'end_rmsd': 0.0, 
        'max_displacement': 0.0,
        'total_path_length': 0.0,
        'smoothness_score': 0.0,
        'avg_step_size': 0.0
    }
    
    # Check endpoint accuracy if reference structures provided
    if start_atoms is not None:
        start_diff = trajectory[0].positions - start_atoms.positions
        metrics['start_rmsd'] = np.sqrt(np.mean(start_diff**2))
        
    if end_atoms is not None and len(trajectory) >= 2:
        end_diff = trajectory[-1].positions - end_atoms.positions  
        metrics['end_rmsd'] = np.sqrt(np.mean(end_diff**2))
    
    # Check path properties
    step_sizes = []
    for i in range(len(trajectory) - 1):
        pos_diff = trajectory[i+1].positions - trajectory[i].positions
        step_displacement = np.sqrt(np.sum(pos_diff**2, axis=1))
        
        max_step = np.max(step_displacement)
        total_step = np.sum(step_displacement)
        
        metrics['max_displacement'] = max(metrics['max_displacement'], max_step)
        metrics['total_path_length'] += total_step
        step_sizes.append(total_step)
    
    if step_sizes:
        metrics['avg_step_size'] = np.mean(step_sizes)
    
    # Calculate smoothness (second derivative measure - lower is smoother)
    if len(trajectory) >= 3:
        second_derivs = []
        for i in range(1, len(trajectory) - 1):
            pos_prev = trajectory[i-1].positions.flatten()
            pos_curr = trajectory[i].positions.flatten()  
            pos_next = trajectory[i+1].positions.flatten()
            
            # Second derivative approximation
            second_deriv = pos_prev - 2*pos_curr + pos_next
            second_derivs.append(np.linalg.norm(second_deriv))
        
        metrics['smoothness_score'] = np.mean(second_derivs)
    
    return metrics


def compare_interpolation_methods(
    start_atoms: Atoms,
    end_atoms: Atoms,
    base_indices: List[int],
    n_frames: int = 20,
    methods: List[str] = None
) -> Dict[str, Dict[str, Any]]:
    """
    Compare different interpolation methods on the same input.
    
    Args:
        start_atoms: Starting structure
        end_atoms: Ending structure
        base_indices: Alignment atoms
        n_frames: Number of frames for comparison
        methods: List of methods to compare (default: ['linear', 'idpp'])
        
    Returns:
        Dictionary with results for each method:
        {
            'linear': {'trajectory': [...], 'quality': {...}, 'time': 0.5},
            'idpp': {'trajectory': [...], 'quality': {...}, 'time': 2.3}
        }
    """
    
    if methods is None:
        methods = ['linear', 'idpp']
    
    results = {}
    
    for method in methods:
        import time
        start_time = time.time()
        
        try:
            # Create interpolated trajectory
            trajectory = create_interpolated_trajectory(
                start_atoms=start_atoms,
                end_atoms=end_atoms,
                n_frames=n_frames,
                method=method,
                base_indices=base_indices
            )
            
            # Assess quality
            quality = interpolation_quality_check(
                trajectory=trajectory,
                start_atoms=start_atoms,
                end_atoms=end_atoms,
                method_name=method
            )
            
            elapsed_time = time.time() - start_time
            
            results[method] = {
                'trajectory': trajectory,
                'quality': quality,
                'time_seconds': elapsed_time,
                'success': True
            }
            
        except Exception as e:
            results[method] = {
                'trajectory': None,
                'quality': None,
                'time_seconds': time.time() - start_time,
                'success': False,
                'error': str(e)
            }
    
    return results