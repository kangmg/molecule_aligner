# Molecule Aligner

A Python package that provides a programmatic API to align and merge multiple ASE trajectories based on a shared base molecule. The trajectories describe ligand binding/unbinding processes and are aligned using specified atom indices of the base molecule.

## Installation

Using uv:
```bash
uv pip install -e .
```

Using pip:
```bash
pip install -e .
```

## Usage

The package provides a main function `align_and_merge_reactions` that aligns and merges multiple reaction trajectories:

```python
from molecule_aligner import align_and_merge_reactions

# Define reaction trajectories
reaction_inputs = [
    {
        'traj_path': 'rxn1.traj',  # Path to trajectory file
        'reverse': False,
        'base_indices': [0, 1, 2, 3, 4]  # Atom indices for alignment
    },
    {
        'traj': [Atoms(), Atoms(), ...],  # Direct ASE Atoms list
        'reverse': True,
        'base_indices': [0, 1, 2, 3, 4]
    }
]

# Align and merge trajectories
merged_traj = align_and_merge_reactions(
    reactions=reaction_inputs,
    output_path='merged.extxyz',
    reference='first'  # or 'reactant'
)
```

## Input Format

Each trajectory is described using a dictionary with the following fields:

- **`traj_path`** (optional): String path to a `.traj` file
- **`traj`** (optional): List of ASE `Atoms` objects
- **`reverse`** (bool): Whether to reverse the trajectory path
- **`base_indices`** (List[int]): Atom indices that correspond to the shared base molecule

**Note:** `traj_path` and `traj` fields are mutually exclusive - only one should be provided.

## Features

- **Trajectory Loading**: Load trajectories from `.traj` files or directly from ASE Atoms lists
- **Flexible Alignment**: Align trajectories using specified atom indices of the base molecule
- **Trajectory Reversal**: Optionally reverse trajectory direction
- **Multiple References**: Choose between 'first' frame or 'reactant' frame as alignment reference
- **Export Support**: Export merged trajectories to `.extxyz` format
- **Error Handling**: Comprehensive validation of input parameters

## API Reference

### `align_and_merge_reactions(reactions, output_path=None, reference='first')`

Main function to align and merge multiple reaction trajectories.

**Parameters:**
- `reactions` (List[Dict]): List of trajectory dictionaries
- `output_path` (str, optional): Path to save merged trajectory as .extxyz
- `reference` (str): Reference frame selection ('first' or 'reactant')

**Returns:**
- `List[Atoms]`: Merged trajectory as list of ASE Atoms objects

**Raises:**
- `ValueError`: For invalid input parameters or missing required fields

## Package Structure

```
molecule_aligner/
├── __init__.py      # Package initialization
├── align.py         # Core alignment and merging logic
├── utils.py         # Alignment utility functions (Kabsch algorithm)
pyproject.toml       # Package metadata and dependencies
```

## Dependencies

- `ase`: Atomic Simulation Environment
- `numpy`: Numerical computations

## Requirements

- Python >= 3.8
- No CLI interface - designed for programmatic use in Jupyter notebooks and Python scripts
- Type hints and comprehensive docstrings included

## Notes

- Uses Kabsch algorithm for optimal molecular alignment
- Supports trajectories of different lengths and structures  
- Aligns all trajectories to the same base orientation before merging
- Prefers `ase.io.read` and `ase.io.write` for file operations