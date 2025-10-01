"""
Unified reaction pathway builder - build_reaction_pathway implementation.
"""

from typing import List, Dict, Union, Optional, Any, Callable
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read, write

from .align import align_and_merge_reactions_enhanced
from .interpolate import create_interpolated_trajectory


def build_reaction_pathway(
    steps: List[Dict[str, Any]],
    base_indices: List[int],
    output_path: Optional[str] = None,
    reference: str = 'first',
    global_settings: Optional[Dict] = None,
    progress_callback: Optional[Callable] = None
) -> List[Atoms]:
    """
    통합 반응 경로 생성기 - 모든 유형의 입력을 하나의 함수로 처리
    
    Args:
        steps: 반응 단계들의 리스트. 각 단계는 다음 중 하나:
            - 'interpolate': 두 구조 사이의 보간
            - 'trajectory': 기존 궤적 사용  
            - 'frame': 단일 프레임 추가
        base_indices: 전역 정렬 기준 원자 인덱스
        output_path: 결과 저장 경로
        reference: 기준 프레임 ('first', 'reactant')
        global_settings: 전역 기본 설정
        progress_callback: 진행상황 콜백 함수(step_index, step_type, progress)
    
    Returns:
        완성된 반응 경로 (List[Atoms])
        
    Example:
        steps = [
            {
                'type': 'interpolate',
                'from': 'start.xyz',
                'to': 'end.xyz', 
                'frames': 20,
                'method': 'idpp'
            },
            {
                'type': 'trajectory',
                'source': 'md.traj',
                'reverse': True,
                'frames': 30
            }
        ]
        
        pathway = build_reaction_pathway(
            steps=steps,
            base_indices=[0,1,2,3,4],
            output_path='complete_pathway.extxyz'
        )
    """
    
    if global_settings is None:
        global_settings = {}
    
    print(f"🚀 Building Reaction Pathway with {len(steps)} steps")
    print(f"🎯 Alignment indices: {len(base_indices)} atoms")
    
    all_reactions = []
    
    for step_idx, step in enumerate(steps):
        step_type = step.get('type', 'unknown')
        
        if progress_callback:
            progress_callback(step_idx, step_type, 0)
        
        print(f"\n📍 Step {step_idx + 1}/{len(steps)}: {step_type}")
        
        try:
            if step_type == 'interpolate':
                reaction = _handle_interpolate_step(step, base_indices, global_settings)
            elif step_type == 'trajectory':
                reaction = _handle_trajectory_step(step, base_indices, global_settings)
            elif step_type == 'frame':
                reaction = _handle_frame_step(step, base_indices, global_settings)
            else:
                raise ValueError(f"Unknown step type: {step_type}")
            
            all_reactions.extend(reaction)
            
            if progress_callback:
                progress_callback(step_idx, step_type, 100)
                
        except Exception as e:
            print(f"❌ Error in step {step_idx + 1}: {e}")
            
            # Try fallback if available
            if step_type == 'interpolate' and step.get('fallback'):
                print(f"🔄 Trying fallback method: {step['fallback']}")
                try:
                    step_fallback = step.copy()
                    step_fallback['method'] = step['fallback']
                    reaction = _handle_interpolate_step(step_fallback, base_indices, global_settings)
                    all_reactions.extend(reaction)
                    print(f"✅ Fallback successful")
                except Exception as e2:
                    print(f"❌ Fallback also failed: {e2}")
                    raise
            else:
                raise
    
    print(f"\n🔧 Assembling final pathway from {len(all_reactions)} segments...")
    
    # Merge all reactions using existing function
    final_trajectory = align_and_merge_reactions_enhanced(
        reactions=all_reactions,
        output_path=output_path,
        reference=reference
    )
    
    print(f"✅ Pathway complete: {len(final_trajectory)} frames")
    
    return final_trajectory


def _handle_interpolate_step(step: Dict, base_indices: List[int], global_settings: Dict) -> List[Dict]:
    """Handle interpolation step with robust type checking."""
    
    # Load source structures with enhanced type safety
    from_source = step['from']
    to_source = step['to']
    
    # Handle 'from' structure
    if isinstance(from_source, str):
        from_atoms = read(from_source)
        if not isinstance(from_atoms, Atoms):
            raise TypeError(f"Loaded 'from' structure is not Atoms object: {type(from_atoms)}")
    elif isinstance(from_source, Atoms):
        from_atoms = from_source
    else:
        raise TypeError(f"'from' must be file path (str) or Atoms object, got: {type(from_source)}")
        
    # Handle 'to' structure  
    if isinstance(to_source, str):
        to_atoms = read(to_source)
        if not isinstance(to_atoms, Atoms):
            raise TypeError(f"Loaded 'to' structure is not Atoms object: {type(to_atoms)}")
    elif isinstance(to_source, Atoms):
        to_atoms = to_source
    else:
        raise TypeError(f"'to' must be file path (str) or Atoms object, got: {type(to_source)}")
    
    # Get parameters with validation
    n_frames = step.get('frames', global_settings.get('default_frames', 10))
    method = step.get('method', global_settings.get('default_method', 'idpp'))
    options = step.get('options', {})
    reverse = step.get('reverse', False)
    
    # Validate n_frames
    if not isinstance(n_frames, int) and n_frames != 'auto':
        raise TypeError(f"n_frames must be int or 'auto', got: {type(n_frames)}")
    
    # Handle 'auto' frame count
    if n_frames == 'auto':
        n_frames = _calculate_auto_frames(from_atoms, to_atoms, base_indices, step)
    
    if n_frames < 2:
        raise ValueError(f"n_frames must be at least 2, got: {n_frames}")
    
    print(f"   🔄 Interpolating {n_frames} frames using {method}")
    print(f"       From: {len(from_atoms)} atoms ({from_atoms.get_chemical_formula()})")
    print(f"       To:   {len(to_atoms)} atoms ({to_atoms.get_chemical_formula()})")
    
    # Create interpolated trajectory
    interpolated_traj = create_interpolated_trajectory(
        start_atoms=from_atoms,
        end_atoms=to_atoms,
        n_frames=n_frames,
        method=method,
        base_indices=base_indices,
        **options
    )
    
    # Validate interpolated trajectory
    if not isinstance(interpolated_traj, list):
        raise TypeError(f"create_interpolated_trajectory returned {type(interpolated_traj)}, expected list")
    
    for i, frame in enumerate(interpolated_traj):
        if not isinstance(frame, Atoms):
            raise TypeError(f"Interpolated frame {i} is not Atoms object: {type(frame)}")
    
    # Apply reverse if requested
    if reverse:
        interpolated_traj = interpolated_traj[::-1]
        print(f"   🔄 Reversed interpolated trajectory")
    
    # Convert to reaction format (all frames from interpolation)
    reactions = []
    for i, frame in enumerate(interpolated_traj):
        reactions.append({
            'traj': [frame],  # Always single-element list for consistency
            'base_indices': base_indices,
            'reverse': False
        })
    
    return reactions


def _handle_trajectory_step(step: Dict, base_indices: List[int], global_settings: Dict) -> List[Dict]:
    """Handle trajectory step with robust type checking."""
    
    source = step['source']
    reverse = step.get('reverse', False)
    frames = step.get('frames', 'all')
    skip = step.get('skip', 1)
    
    # Load trajectory with robust type checking
    if isinstance(source, str):
        trajectory = read(source, ':')
        # Ensure trajectory is always a list
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
    elif isinstance(source, list):
        # Verify it's a list of Atoms objects
        if not all(isinstance(item, Atoms) for item in source):
            raise ValueError("trajectory source list must contain only Atoms objects")
        trajectory = source
    else:
        raise ValueError("Trajectory source must be file path (str) or list of Atoms objects")
    
    print(f"   📂 Loaded trajectory: {len(trajectory)} frames")
    
    # Handle frame selection
    if frames != 'all':
        if isinstance(frames, int):
            # Take first N frames
            trajectory = trajectory[:frames]
        elif isinstance(frames, list) and len(frames) == 2:
            # Take frame range [start, end]
            start, end = frames
            trajectory = trajectory[start:end]
        else:
            raise ValueError("frames must be 'all', int, or [start, end]")
    
    # Handle frame skipping
    if skip > 1:
        trajectory = trajectory[::skip]
        print(f"   ⏭️  Skipped frames, now {len(trajectory)} frames")
    
    # Apply reverse
    if reverse:
        trajectory = trajectory[::-1]
        print(f"   🔄 Reversed trajectory")
    
    # Convert to reaction format with type safety
    reactions = []
    for i, frame in enumerate(trajectory):
        if not isinstance(frame, Atoms):
            raise TypeError(f"Frame {i} is not an Atoms object: {type(frame)}")
        
        reactions.append({
            'traj': [frame],  # Always wrap single frame in list for consistency
            'base_indices': base_indices,
            'reverse': False
        })
    
    return reactions


def _handle_frame_step(step: Dict, base_indices: List[int], global_settings: Dict) -> List[Dict]:
    """Handle single frame step."""
    
    source = step['source']
    repeat = step.get('repeat', 1)
    
    # Load frame
    if isinstance(source, str):
        frame = read(source)
    else:
        frame = source
    
    print(f"   📄 Added frame, repeat {repeat} times")
    
    # Create reactions with repetition
    reactions = []
    for _ in range(repeat):
        reactions.append({
            'traj': [frame.copy()],
            'base_indices': base_indices,
            'reverse': False
        })
    
    return reactions


def _calculate_auto_frames(from_atoms: Atoms, to_atoms: Atoms, base_indices: List[int], step: Dict) -> int:
    """Calculate automatic frame count based on structural difference."""
    
    # Get alignment positions
    pos1 = from_atoms.positions[base_indices]
    pos2 = to_atoms.positions[base_indices]
    
    # Calculate RMSD
    centroid1 = np.mean(pos1, axis=0)
    centroid2 = np.mean(pos2, axis=0)
    
    pos2_aligned = pos2 - centroid2 + centroid1
    rmsd = np.sqrt(np.mean((pos1 - pos2_aligned)**2))
    
    # Get frames per angstrom setting
    frames_per_angstrom = step.get('frames_per_angstrom', 5)
    
    # Calculate frame count
    calculated_frames = max(5, int(rmsd * frames_per_angstrom))
    
    print(f"   📐 Auto frames: RMSD={rmsd:.2f}Å → {calculated_frames} frames")
    
    return calculated_frames


# Convenience functions for common patterns
def create_interpolation_chain(structures: List[Union[str, Atoms]], 
                              base_indices: List[int],
                              frames_per_segment: Union[int, List[int]] = 15,
                              method: str = 'idpp',
                              **kwargs) -> List[Atoms]:
    """
    Create a chain of interpolations: A → B → C → D ...
    
    Args:
        structures: List of structures (file paths or Atoms objects)
        base_indices: Alignment indices
        frames_per_segment: Frames for each interpolation (int or list)
        method: Interpolation method
        **kwargs: Additional arguments for build_reaction_pathway
    
    Returns:
        Complete interpolated chain
    """
    
    if isinstance(frames_per_segment, int):
        frames_list = [frames_per_segment] * (len(structures) - 1)
    else:
        frames_list = frames_per_segment
    
    if len(frames_list) != len(structures) - 1:
        raise ValueError("frames_per_segment length must match number of segments")
    
    steps = []
    for i in range(len(structures) - 1):
        steps.append({
            'type': 'interpolate',
            'from': structures[i],
            'to': structures[i + 1],
            'frames': frames_list[i],
            'method': method
        })
    
    return build_reaction_pathway(
        steps=steps,
        base_indices=base_indices,
        **kwargs
    )


def create_cyclic_pathway(structures: List[Union[str, Atoms]],
                         base_indices: List[int], 
                         frames_per_segment: int = 15,
                         method: str = 'idpp',
                         **kwargs) -> List[Atoms]:
    """
    Create cyclic pathway: A → B → C → A
    
    Args:
        structures: List of structures (will be made cyclic)
        base_indices: Alignment indices
        frames_per_segment: Frames for each segment
        method: Interpolation method
        **kwargs: Additional arguments
    
    Returns:
        Cyclic pathway
    """
    
    # Make cyclic by adding first structure at the end
    cyclic_structures = structures + [structures[0]]
    
    return create_interpolation_chain(
        structures=cyclic_structures,
        base_indices=base_indices,
        frames_per_segment=frames_per_segment,
        method=method,
        **kwargs
    )