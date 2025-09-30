# qwen.md

## Task

Write a Python package that provides a programmatic API to align and merge multiple ASE trajectories based on a shared base molecule. The trajectories describe ligand binding/unbinding processes, and should be aligned using specified atom indices of the base molecule.

## Usage Context

The API is intended to be used in Jupyter notebooks and Python scripts, not from the command line.

## Input Format

Each trajectory is described using a dictionary with the following fields:

- '''traj_path''' (optional): string path to a ''.traj'' file.
- '''traj''' (optional): list of ASE '''Atoms''' objects.
- '''reverse''' (bool): whether to reverse the trajectory path.
- '''base_indices''' (List[int]): atom indices that correspond to the shared base molecule.

Example input:

'''
reaction_inputs = [
    {
        'traj_path': 'rxn1.traj',
        'reverse': False,
        'base_indices': [0, 1, 2, 3, 4]
    },
    {
        'traj': [Atoms(), Atoms(), ...],
        'reverse': True,
        'base_indices': [0, 1, 2, 3, 4]
    }
]
'''

## Function Signature

Your main API should look like this:

'''
from molecule_aligner import align_and_merge_reactions

merged_traj = align_and_merge_reactions(
    reactions=reaction_inputs,
    output_path='merged.extxyz',
    reference='first'  # or 'reactant'
)
'''

## Requirements

- Accept a list of dicts as input, each representing a reaction trajectory.
- Each reaction trajectory must be aligned using the '''base_indices''' provided.
- Alignment must superimpose the '''base_indices''' onto the reference frame.
- The function should:
  - Optionally reverse the trajectory if '''reverse=True'''.
  - Align all trajectories to the same base orientation.
  - Merge all aligned frames into a single trajectory.
  - Export the final result to ''.extxyz''.
  - Return the merged ASE '''Trajectory'''.

## Notes

- Use '''ase.geometry''' or RMSD-based alignment with NumPy.
- Prefer '''ase.io.read''' and '''ase.io.write''' for file operations.
- Do not assume all trajectories are the same length or structure.
- The '''traj''' and '''traj_path''' fields are mutually exclusive; either one may be provided.

## Package Structure

Generate a Python package with the following structure:

'''
molecule_aligner/
├── __init__.py
├── align.py        # core alignment and merging logic
├── utils.py        # optional helper functions
pyproject.toml      # package metadata
'''

## Implementation Constraints

- Use '''pyproject.toml''' for build system and dependencies ('''ase''', '''numpy''').
- Do not create a CLI interface.
- All functions must be importable and callable from Python code.
- Add type hints and docstrings.
- Use '''overlay.py''' as optional reference — re-implement or borrow logic as appropriate.

