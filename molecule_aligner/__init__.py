"""
Molecule Aligner Package

A Python package to align and merge multiple ASE trajectories based on a shared base molecule.
"""

from .align import align_and_merge_reactions

__version__ = "0.1.0"
__all__ = ["align_and_merge_reactions"]