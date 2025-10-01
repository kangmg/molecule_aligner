# 🧬 Molecule Aligner

A powerful Python package for creating complex molecular reaction pathways by seamlessly combining trajectory alignment, interpolation, and frame manipulation in a single, intuitive API.

## 🚀 Overview

Molecule Aligner enables computational chemists to build sophisticated molecular reaction mechanisms from diverse sources - existing MD trajectories, single molecular structures, and interpolated transition paths. Whether you're studying catalytic cycles, protein conformational changes, or drug binding mechanisms, Molecule Aligner provides the tools to create smooth, chemically accurate reaction pathways.

## ✨ Key Features

- **🎯 Unified API**: Single function handles all pathway construction needs
- **🧪 Mixed Workflows**: Seamlessly combine trajectories, interpolations, and static frames
- **⚗️ Chemical Accuracy**: IDPP interpolation for realistic transition paths
- **🔄 Flexible Input**: Support for files, ASE objects, and various formats
- **🛡️ Robust Operation**: Automatic fallbacks and comprehensive error handling
- **📊 Quality Control**: Built-in assessment tools and method comparison

## 🔧 Installation

### Using uv (recommended)

```bash
git clone <repository-url>
cd molecule_aligner
uv venv && source .venv/bin/activate
uv pip install -e .
```

### Using pip

```bash
pip install -e .
```

## 🎯 Quick Start

### Basic Reaction Pathway

```python
from molecule_aligner import build_reaction_pathway

# Create a complete reaction pathway with mixed inputs
pathway = build_reaction_pathway(
    steps=[
        # Interpolate initial approach
        {
            'type': 'interpolate',
            'from': 'reactant.xyz',
            'to': 'intermediate.xyz',
            'frames': 15,
            'method': 'idpp'
        },
        
        # Use existing MD trajectory  
        {
            'type': 'trajectory',
            'source': 'reaction_dynamics.traj',
            'frames': 'all',
            'reverse': False
        },
        
        # Final product formation
        {
            'type': 'interpolate', 
            'from': 'intermediate.xyz',
            'to': 'product.xyz',
            'frames': 12,
            'method': 'idpp'
        }
    ],
    base_indices=[0, 1, 2, 3, 4, 5],  # Atoms for alignment
    output_path='complete_reaction.extxyz'
)

print(f"Generated pathway: {len(pathway)} frames")
```

### Catalytic Cycle Example

```python
# Complex catalytic mechanism
catalytic_cycle = build_reaction_pathway(
    steps=[
        # Substrate binding
        {'type': 'interpolate', 'from': 'catalyst.xyz', 'to': 'bound.xyz', 'frames': 20, 'method': 'idpp'},
        
        # Reaction dynamics (MD trajectory)  
        {'type': 'trajectory', 'source': 'reaction_md.traj', 'frames': [100, 300], 'skip': 2},
        
        # Hold at transition state
        {'type': 'frame', 'source': 'transition_state.xyz', 'repeat': 5},
        
        # Product release
        {'type': 'interpolate', 'from': 'product_bound.xyz', 'to': 'catalyst.xyz', 'frames': 15, 'method': 'idpp'}
    ],
    base_indices=list(range(20)),
    output_path='catalytic_cycle.extxyz'
)
```

## 🎮 Step Types

Molecule Aligner supports three fundamental step types for maximum flexibility:

### `interpolate` - Smooth Transitions
Create chemically accurate transitions between molecular structures using IDPP or linear methods.

### `trajectory` - Existing Simulations  
Incorporate MD trajectories, conformational sampling, or any multi-frame molecular data.

### `frame` - Static States
Add checkpoints, static periods, or important intermediate structures.

> 📚 **For complete step reference and advanced examples, see [GUIDE.md](GUIDE.md)**

## 🧪 Real-World Applications

- **🔬 Reaction Mechanisms**: Complete multi-step organic and organometallic reactions
- **💊 Drug Discovery**: Detailed protein-ligand binding pathways  
- **🧬 Protein Dynamics**: Conformational changes and folding pathways
- **⚡ Catalysis**: Full catalytic cycles from substrate binding to product release
- **🔄 Reversible Processes**: Bidirectional reactions and equilibrium dynamics

## 🛠️ Interpolation Methods

### IDPP (Image Dependent Pair Potential)
- **Accuracy**: Chemically realistic transition paths
- **Use case**: High-quality reaction mechanisms
- **Performance**: Slower but more accurate

### Linear Interpolation  
- **Speed**: Fast geometric interpolation
- **Use case**: Prototyping and simple transitions
- **Performance**: Very fast but less chemically accurate

### Automatic Fallback
Robust operation with automatic fallback from IDPP to linear if needed.

## 📁 Input/Output Formats

### Supported Inputs
- ASE trajectory files (`.traj`)
- XYZ files (`.xyz`) 
- PDB files (`.pdb`)
- Any ASE-supported format
- Direct `Atoms` objects
- Lists of `Atoms` objects

### Output Formats
- Extended XYZ with metadata (`.extxyz`)
- ASE trajectory format (`.traj`)
- Automatic format detection

## 📊 Visualization

```bash
# ASE GUI (recommended)
ase gui your_pathway.extxyz

# VMD
vmd your_pathway.extxyz  

# PyMOL (convert first)
ase convert your_pathway.extxyz your_pathway.pdb
```

## 🎓 Documentation

- **[GUIDE.md](GUIDE.md)** - Complete step reference and advanced examples
- **[ADVANCED_USAGE_GUIDE.md](ADVANCED_USAGE_GUIDE.md)** - Real-world applications and best practices
- **[UNIFIED_API_PROPOSAL.md](UNIFIED_API_PROPOSAL.md)** - Technical design and future features

## 🚀 Performance Tips

- Use `'linear'` method for prototyping (1000x faster)
- Use `'idpp'` for final high-quality results  
- Start with fewer frames, increase for final production
- Use `skip` parameter for large trajectories
- Test with small systems first

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch  
3. Add tests for new functionality
4. Submit a pull request

## 📄 License

MIT License - see LICENSE file for details.

## 📚 Citation

```bibtex
@software{molecule_aligner,
  title={Molecule Aligner: Unified Molecular Pathway Construction},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/molecule_aligner}
}
```

## 🔄 Version History

### v0.2.0 (Current)
- ✨ **NEW**: `build_reaction_pathway()` unified API
- 🧪 IDPP interpolation for chemical accuracy
- 🔄 Mixed workflow support (trajectories + interpolations + frames)
- 🛡️ Enhanced error handling and fallbacks
- 📚 Comprehensive documentation and examples

### v0.1.0
- Basic trajectory alignment and merging functionality

---

**Ready to build your molecular pathways? Check out [GUIDE.md](GUIDE.md) for detailed examples!** 🚀