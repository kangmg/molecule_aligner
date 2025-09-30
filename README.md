# Molecule Aligner

A Python package that provides a programmatic API to align and merge multiple ASE trajectories based on a shared base molecule. Enhanced with **single frame interpolation capabilities** using ASE's IDPP method for generating chemically accurate reaction paths.

## 🆕 Version 0.2.0 - New Features!

- ✨ **Single Frame Interpolation**: Generate reaction paths from just two structures
- 🧪 **IDPP Method**: Uses ASE's Image Dependent Pair Potential for chemically realistic paths  
- 🔀 **Mixed Workflows**: Combine existing trajectories with interpolated segments
- 📊 **Quality Assessment**: Built-in metrics for path analysis
- 🎯 **Convenient API**: Simple functions for common use cases
- 🔄 **Backward Compatible**: All existing code continues to work

## Installation

Using uv:
```bash
uv pip install -e .
```

Using pip:
```bash
pip install -e .
```

## Quick Start Examples

### 🚀 Simple Reaction Path Generation (NEW!)
```python
from molecule_aligner import create_reaction_path

# Generate chemically accurate reaction path between two structures
trajectory = create_reaction_path(
    start_structure='reactant.xyz',
    end_structure='product.xyz', 
    base_indices=[0, 1, 2, 3],    # Atoms to keep aligned
    n_frames=20,                  # Number of intermediate frames
    method='idpp',                # Use IDPP for chemical accuracy
    output_path='reaction_path.extxyz'
)

print(f"Generated {len(trajectory)} frame reaction path!")
```

### 🔬 Method Comparison
```python
from molecule_aligner.interpolate import compare_interpolation_methods

# Compare different interpolation methods
comparison = compare_interpolation_methods(
    start_atoms=reactant,
    end_atoms=product,
    base_indices=[0, 1, 2],
    methods=['linear', 'idpp']
)

for method, results in comparison.items():
    quality = results['quality']
    print(f"{method}: smoothness={quality['smoothness_score']:.3f}")
```

### 🔄 Mixed Trajectory + Interpolation Workflow
```python
from molecule_aligner import align_and_merge_reactions_enhanced

reactions = [
    {
        # Existing multi-frame trajectory
        'traj_path': 'binding.traj',
        'reverse': False,
        'base_indices': [0, 1, 2, 3, 4, 5]
    },
    {
        # NEW: Single frame interpolation  
        'start_frame_path': 'intermediate.xyz',
        'end_frame_path': 'product.xyz',
        'n_frames': 15,
        'interpolation_method': 'idpp',
        'interpolation_options': {
            'fmax': 0.05,              # Tight convergence
            'optimizer': 'LBFGS'       # High-quality optimizer
        },
        'base_indices': [0, 1, 2, 3, 4, 5]
    }
]

complete_trajectory = align_and_merge_reactions_enhanced(reactions)
```

## Input Formats

### Traditional Trajectory Format (Existing)
```python
{
    'traj_path': 'trajectory.traj',    # OR 'traj': [Atoms(), ...]
    'reverse': False,
    'base_indices': [0, 1, 2, 3, 4]
}
```

### NEW: Single Frame Interpolation Format
```python
{
    'start_frame': Atoms(),            # OR 'start_frame_path': 'start.xyz'
    'end_frame': Atoms(),              # OR 'end_frame_path': 'end.xyz'
    'n_frames': 20,                    # Number of interpolated frames
    'interpolation_method': 'idpp',    # 'linear' or 'idpp' 
    'interpolation_options': {         # Method-specific options
        'fmax': 0.1,                   # IDPP convergence criterion
        'steps': 100,                  # IDPP optimization steps
        'optimizer': 'LBFGS',          # 'MDMin', 'LBFGS', or 'FIRE'
        'mic': False                   # Minimum image convention
    },
    'base_indices': [0, 1, 2, 3, 4],
    'reverse': False
}
```

## API Reference

### Enhanced Functions (v0.2.0)

#### `create_reaction_path(start_structure, end_structure, base_indices, ...)`
Simple interface for generating reaction paths between two structures.

#### `align_and_merge_reactions_enhanced(reactions, ...)`
Enhanced version supporting both trajectories and interpolation.

#### Quality Assessment
```python
from molecule_aligner.interpolate import interpolation_quality_check

quality = interpolation_quality_check(trajectory, start_atoms, end_atoms)
print(f"Path smoothness: {quality['smoothness_score']}")
print(f"Total path length: {quality['total_path_length']} Å")
```

### Original Functions (Backward Compatible)

#### `align_and_merge_reactions(reactions, ...)`
Original function - works exactly as before in v0.1.0.

## Interpolation Methods

### Linear Interpolation
- **Speed**: Very fast (~1000 frames/second)
- **Use Case**: Simple transitions, quick previews
- **Quality**: Basic, may create unrealistic intermediate structures

### IDPP (Image Dependent Pair Potential)  
- **Speed**: Moderate (~10-100 frames/second)
- **Use Case**: Chemically accurate reaction paths
- **Quality**: High - considers interatomic distances and avoids unphysical configurations
- **Optimizers**: MDMin (default), LBFGS, FIRE

## Quality Metrics

The package provides comprehensive quality assessment:

```python
quality_metrics = {
    'n_frames': 20,
    'start_rmsd': 0.0,              # Accuracy of start point
    'end_rmsd': 0.0,                # Accuracy of end point  
    'max_displacement': 0.5,        # Largest single step
    'total_path_length': 15.2,      # Total path length (Å)
    'smoothness_score': 0.8,        # Path smoothness (lower = smoother)
    'avg_step_size': 0.76           # Average step size (Å)
}
```

## Performance Characteristics

| Method | Speed | Memory | Chemical Accuracy | Use Case |
|--------|-------|---------|-------------------|----------|
| Linear | ⚡⚡⚡ | ✅ | ⭐ | Quick previews |
| IDPP | ⚡⚡ | ✅ | ⭐⭐⭐ | Reaction paths |

## Dependencies

- **ase** >= 3.22.0 (includes IDPP implementation)
- **numpy** >= 1.20.0

## Migration from v0.1.0

**No changes required!** All existing code works exactly as before:

```python
# This still works identically
from molecule_aligner import align_and_merge_reactions

result = align_and_merge_reactions(reactions)  # Same as always
```

**To use new features**, import enhanced functions:

```python  
# New interpolation capabilities
from molecule_aligner import create_reaction_path, align_and_merge_reactions_enhanced
```

## Examples and Use Cases

### Protein-Ligand Binding
```python
# Model ligand approach to binding site
binding_path = create_reaction_path(
    start_structure='ligand_unbound.pdb',
    end_structure='ligand_bound.pdb',
    base_indices=[0, 1, 2, 3, 4],  # Protein backbone atoms
    method='idpp',
    n_frames=30
)
```

### Conformational Changes
```python
# Model protein conformational transition
conformational_path = create_reaction_path(
    start_structure=closed_conformation,
    end_structure=open_conformation, 
    base_indices=active_site_atoms,
    method='idpp',
    fmax=0.02  # High accuracy
)
```

### Multi-Step Reaction Pathways
```python
# Complex reaction with multiple intermediates
reactions = [
    {'start_frame': reactant, 'end_frame': intermediate1, 'n_frames': 10, 'interpolation_method': 'idpp'},
    {'start_frame': intermediate1, 'end_frame': intermediate2, 'n_frames': 8, 'interpolation_method': 'idpp'}, 
    {'start_frame': intermediate2, 'end_frame': product, 'n_frames': 12, 'interpolation_method': 'idpp'}
]
```

## Notes

- **IDPP Method**: Leverages ASE's proven implementation for production-quality results
- **Memory Efficient**: Processes large molecules and long trajectories efficiently  
- **Error Handling**: Graceful fallbacks and comprehensive validation
- **File Format Support**: Reads various formats (.xyz, .traj, .pdb, etc.) via ASE
- **Visualization Ready**: Output compatible with VMD, Ovito, ASE GUI