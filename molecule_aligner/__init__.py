"""
Molecule Aligner Package

A Python package to align and merge multiple ASE trajectories based on a shared base molecule.
Enhanced with single frame interpolation capabilities using ASE IDPP method.
"""

from .align import (
    align_and_merge_reactions,           # Original function (backward compatibility)
    align_and_merge_reactions_enhanced,  # Enhanced function with interpolation
    create_reaction_path                 # Convenience function for interpolation
)

__version__ = "0.2.0"
__all__ = [
    "align_and_merge_reactions",
    "align_and_merge_reactions_enhanced", 
    "create_reaction_path"
]