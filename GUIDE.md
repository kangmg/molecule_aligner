# 📚 Molecule Aligner Complete Guide

## 🚀 Overview

Molecule Aligner provides a unified API for creating complex molecular reaction pathways by combining trajectory alignment, interpolation, and frame manipulation in a single function call.

## 🎯 Core Function: `build_reaction_pathway()`

```python
from molecule_aligner import build_reaction_pathway

pathway = build_reaction_pathway(
    steps: List[Dict],              # List of step definitions
    base_indices: List[int],        # Atoms for alignment (e.g., [0,1,2,3,4])
    output_path: str = None,        # Output file path (optional)
    reference: str = 'first',       # Reference frame ('first', 'reactant')
    global_settings: Dict = None,   # Global default settings
    progress_callback: Callable = None  # Progress monitoring function
) -> List[Atoms]
```

## 📋 Step Types Reference

### 🔹 1. `interpolate` - Molecular Interpolation

Creates smooth transitions between two molecular structures using linear or IDPP methods.

#### **Parameters**

| Parameter | Type | Required | Description | Default |
|-----------|------|----------|-------------|---------|
| `type` | `str` | ✅ Yes | Must be `'interpolate'` | - |
| `from` | `str` or `Atoms` | ✅ Yes | Start structure (file path or Atoms object) | - |
| `to` | `str` or `Atoms` | ✅ Yes | End structure (file path or Atoms object) | - |
| `frames` | `int` or `'auto'` | ✅ Yes | Number of interpolation frames (≥2) | - |
| `method` | `str` | No | Interpolation method: `'linear'` or `'idpp'` | `'idpp'` |
| `options` | `Dict` | No | Method-specific options (see below) | `{}` |
| `reverse` | `bool` | No | Reverse the interpolated trajectory | `False` |
| `fallback` | `str` | No | Fallback method if primary fails | `None` |

#### **Method Options**

**For `'idpp'` method:**
```python
'options': {
    'fmax': 0.1,           # Convergence criterion (force tolerance)
    'steps': 100,          # Maximum optimization steps
    'optimizer': 'LBFGS',  # Optimizer: 'LBFGS', 'FIRE', 'MDMin'
    'mic': False           # Minimum image convention for PBC
}
```

**For `'linear'` method:**
```python
'options': {
    'mic': False           # Minimum image convention for PBC
}
```

#### **Examples**

**Basic IDPP interpolation:**
```python
{
    'type': 'interpolate',
    'from': 'reactant.xyz',
    'to': 'product.xyz', 
    'frames': 20,
    'method': 'idpp'
}
```

**High-precision interpolation:**
```python
{
    'type': 'interpolate',
    'from': start_atoms,
    'to': end_atoms,
    'frames': 30,
    'method': 'idpp',
    'options': {
        'fmax': 0.05,      # Higher precision
        'steps': 200,      # More optimization steps
        'optimizer': 'LBFGS'
    }
}
```

**Fast linear interpolation:**
```python
{
    'type': 'interpolate',
    'from': 'intermediate1.xyz',
    'to': 'intermediate2.xyz',
    'frames': 8,
    'method': 'linear'     # Fast but less accurate
}
```

**Auto frame count:**
```python
{
    'type': 'interpolate',
    'from': mol1,
    'to': mol2,
    'frames': 'auto',      # Automatically determine based on RMSD
    'frames_per_angstrom': 5,  # 5 frames per Å of movement
    'method': 'idpp'
}
```

**With fallback mechanism:**
```python
{
    'type': 'interpolate',
    'from': 'complex_start.xyz',
    'to': 'complex_end.xyz',
    'frames': 15,
    'method': 'idpp',
    'fallback': 'linear',  # Use linear if IDPP fails
    'options': {'fmax': 0.1}
}
```

**Reversed interpolation:**
```python
{
    'type': 'interpolate',
    'from': 'product.xyz',
    'to': 'reactant.xyz',
    'frames': 12,
    'method': 'idpp',
    'reverse': True        # Reverse the trajectory order
}
```

---

### 🔹 2. `trajectory` - Trajectory Processing

Uses existing molecular dynamics trajectories or frame sequences.

#### **Parameters**

| Parameter | Type | Required | Description | Default |
|-----------|------|----------|-------------|---------|
| `type` | `str` | ✅ Yes | Must be `'trajectory'` | - |
| `source` | `str` or `List[Atoms]` | ✅ Yes | Trajectory file path or list of Atoms | - |
| `frames` | `str`, `int`, or `List[int]` | ✅ Yes | Frame selection (see options below) | - |
| `reverse` | `bool` | No | Reverse trajectory order | `False` |
| `skip` | `int` | No | Frame skipping (1=all, 2=every other, etc.) | `1` |

#### **Frame Selection Options**

| `frames` Value | Description | Example |
|----------------|-------------|---------|
| `'all'` | Use all frames from trajectory | `'frames': 'all'` |
| `int` | Use first N frames | `'frames': 50` |
| `[start, end]` | Use frame range | `'frames': [100, 200]` |

#### **Examples**

**Complete MD trajectory:**
```python
{
    'type': 'trajectory',
    'source': 'md_simulation.traj',
    'frames': 'all',
    'reverse': False
}
```

**List of Atoms objects:**
```python
trajectory_list = [atoms1, atoms2, atoms3, atoms4]
{
    'type': 'trajectory',
    'source': trajectory_list,    # Direct List[Atoms]
    'frames': 'all',
    'reverse': False
}
```

**Limited frames:**
```python
{
    'type': 'trajectory',
    'source': 'long_simulation.traj',
    'frames': 100,               # Only first 100 frames
    'skip': 2,                   # Every 2nd frame (50 frames total)
    'reverse': False
}
```

**Frame range selection:**
```python
{
    'type': 'trajectory',
    'source': 'equilibration.traj',
    'frames': [500, 1000],       # Frames 500-1000 only
    'skip': 5,                   # Every 5th frame
    'reverse': False
}
```

**Reversed trajectory (e.g., for reverse reaction):**
```python
{
    'type': 'trajectory',
    'source': 'forward_reaction.traj',
    'frames': 'all',
    'reverse': True              # Reverse order for backward reaction
}
```

**Sampling from long trajectory:**
```python
{
    'type': 'trajectory',
    'source': '/path/to/huge_trajectory.traj',
    'frames': [1000, 5000],      # Take middle section
    'skip': 10,                  # Sample every 10th frame
    'reverse': False
}
```

---

### 🔹 3. `frame` - Single Frame Insertion

Adds individual molecular structures, useful for checkpoints or static states.

#### **Parameters**

| Parameter | Type | Required | Description | Default |
|-----------|------|----------|-------------|---------|
| `type` | `str` | ✅ Yes | Must be `'frame'` | - |
| `source` | `str` or `Atoms` | ✅ Yes | Single structure (file path or Atoms object) | - |
| `repeat` | `int` | No | Number of times to repeat the frame | `1` |

#### **Examples**

**Single checkpoint:**
```python
{
    'type': 'frame',
    'source': 'important_intermediate.xyz',
    'repeat': 1              # Add once
}
```

**Hold at structure (static period):**
```python
{
    'type': 'frame',
    'source': transition_state_atoms,
    'repeat': 10             # Hold for 10 frames
}
```

**Equilibration pause:**
```python
{
    'type': 'frame',
    'source': 'equilibrated_structure.xyz',
    'repeat': 5              # 5-frame pause for equilibration
}
```

**Atoms object insertion:**
```python
critical_geometry = Atoms('H2O', positions=[[0,0,0], [1,0,0], [0,1,0]])
{
    'type': 'frame',
    'source': critical_geometry,  # Direct Atoms object
    'repeat': 3
}
```

---

## 🎮 Complete Workflow Examples

### 🧪 Example 1: Basic Reaction Pathway

```python
from molecule_aligner import build_reaction_pathway

# Simple A → B → C reaction
pathway = build_reaction_pathway(
    steps=[
        # Step 1: Reactant to intermediate
        {
            'type': 'interpolate',
            'from': 'reactant_A.xyz',
            'to': 'intermediate_B.xyz',
            'frames': 15,
            'method': 'idpp'
        },
        # Step 2: Intermediate to product  
        {
            'type': 'interpolate',
            'from': 'intermediate_B.xyz',
            'to': 'product_C.xyz',
            'frames': 12,
            'method': 'idpp'
        }
    ],
    base_indices=[0, 1, 2, 3, 4],
    output_path='A_to_C_pathway.extxyz'
)

print(f"Generated pathway: {len(pathway)} frames")
```

### 🧪 Example 2: Complex Catalytic Cycle

```python
# Comprehensive catalytic cycle with multiple methods
catalytic_cycle = build_reaction_pathway(
    steps=[
        # 1. Substrate approach (IDPP interpolation)
        {
            'type': 'interpolate',
            'from': 'free_catalyst.xyz',
            'to': 'substrate_bound.xyz',
            'frames': 20,
            'method': 'idpp',
            'options': {'fmax': 0.05}  # High precision
        },
        
        # 2. Brief stabilization (static frame)
        {
            'type': 'frame',
            'source': 'substrate_bound.xyz',
            'repeat': 3
        },
        
        # 3. Reaction dynamics (MD trajectory)
        {
            'type': 'trajectory',
            'source': 'reaction_md.traj',
            'frames': [100, 300],      # Use frames 100-300
            'skip': 2,                 # Every 2nd frame
            'reverse': False
        },
        
        # 4. Product formation (linear interpolation)
        {
            'type': 'interpolate',
            'from': 'reaction_intermediate.xyz',
            'to': 'product_bound.xyz',
            'frames': 10,
            'method': 'linear'         # Faster method
        },
        
        # 5. Product release (IDPP)
        {
            'type': 'interpolate',
            'from': 'product_bound.xyz',
            'to': 'regenerated_catalyst.xyz',
            'frames': 15,
            'method': 'idpp'
        }
    ],
    base_indices=list(range(20)),      # First 20 atoms for alignment
    output_path='complete_catalytic_cycle.extxyz'
)
```

### 🧪 Example 3: Protein Conformational Change

```python
# Protein folding pathway with mixed approaches
protein_folding = build_reaction_pathway(
    steps=[
        # 1. Initial unfolding (many frames for detail)
        {
            'type': 'interpolate',
            'from': 'native_protein.pdb',
            'to': 'partially_unfolded.pdb',
            'frames': 50,
            'method': 'idpp',
            'options': {
                'fmax': 0.2,           # Relaxed for flexible protein
                'optimizer': 'FIRE',   # Good for proteins
                'steps': 150
            }
        },
        
        # 2. Unfolded state dynamics (MD simulation)
        {
            'type': 'trajectory',
            'source': 'unfolded_dynamics.traj',
            'frames': 'all',
            'skip': 5,                 # Sample every 5th frame
            'reverse': False
        },
        
        # 3. Refolding initiation (key checkpoint)
        {
            'type': 'frame',
            'source': 'refolding_nucleus.pdb',
            'repeat': 5                # Emphasize important state
        },
        
        # 4. Final refolding (auto frame count)
        {
            'type': 'interpolate',
            'from': 'refolding_nucleus.pdb',
            'to': 'native_protein.pdb',
            'frames': 'auto',
            'frames_per_angstrom': 3,  # 3 frames per Å movement
            'method': 'idpp'
        }
    ],
    base_indices=list(range(30, 50)),  # Use conserved core residues
    output_path='protein_folding_pathway.extxyz'
)
```

### 🧪 Example 4: Drug Binding Process

```python
# Detailed drug-target binding mechanism
drug_binding = build_reaction_pathway(
    steps=[
        # 1. Long-range approach (linear - simple diffusion)
        {
            'type': 'interpolate',
            'from': 'drug_distant.xyz',
            'to': 'drug_surface.xyz', 
            'frames': 8,
            'method': 'linear'
        },
        
        # 2. Surface exploration (experimental MD data)
        {
            'type': 'trajectory',
            'source': exploration_frames,  # List[Atoms] from simulation
            'frames': 'all',
            'reverse': False
        },
        
        # 3. Induced fit binding (high-detail IDPP)
        {
            'type': 'interpolate',
            'from': 'drug_surface.xyz',
            'to': 'drug_bound_open.xyz',
            'frames': 35,
            'method': 'idpp',
            'options': {
                'fmax': 0.02,          # Very high precision
                'steps': 300           # Many optimization steps
            },
            'fallback': 'linear'       # Safety fallback
        },
        
        # 4. Binding site closure (reverse previous opening)
        {
            'type': 'trajectory',
            'source': 'binding_site_closure.traj',
            'frames': [50, 150],
            'reverse': False
        },
        
        # 5. Final optimization (brief)
        {
            'type': 'interpolate',
            'from': 'drug_bound_open.xyz',
            'to': 'drug_bound_closed.xyz',
            'frames': 6,
            'method': 'linear'
        }
    ],
    base_indices=list(range(100, 150)), # Active site residues only
    output_path='drug_binding_mechanism.extxyz',
    
    # Global settings for consistency
    global_settings={
        'default_method': 'idpp',
        'default_frames': 15
    }
)
```

### 🧪 Example 5: Reaction with Reverse Path

```python
# Forward and reverse reaction in one pathway
bidirectional_reaction = build_reaction_pathway(
    steps=[
        # Forward reaction: A → B
        {
            'type': 'interpolate',
            'from': 'reactant_A.xyz',
            'to': 'product_B.xyz',
            'frames': 25,
            'method': 'idpp'
        },
        
        # Equilibrium dynamics at product
        {
            'type': 'trajectory',
            'source': 'product_equilibrium.traj',
            'frames': [0, 50],         # First 50 frames
            'skip': 2
        },
        
        # Reverse reaction: B → A (using reverse)
        {
            'type': 'interpolate',
            'from': 'product_B.xyz',
            'to': 'reactant_A.xyz',
            'frames': 20,
            'method': 'idpp',
            'reverse': True            # This will be reversed again
        }
    ],
    base_indices=[0, 1, 2, 3, 4, 5],
    output_path='bidirectional_reaction.extxyz'
)
```

---

## 🔧 Advanced Features

### 📊 Progress Monitoring

```python
def my_progress_callback(step_index, step_type, progress_percent):
    print(f"Processing step {step_index + 1} ({step_type}): {progress_percent}%")

pathway = build_reaction_pathway(
    steps=my_steps,
    base_indices=[0, 1, 2, 3],
    progress_callback=my_progress_callback  # Monitor progress
)
```

### 🛡️ Error Handling and Fallbacks

```python
# Robust interpolation with multiple fallbacks
robust_steps = [
    {
        'type': 'interpolate',
        'from': complex_start,
        'to': complex_end,
        'frames': 20,
        'method': 'idpp',
        'fallback': 'linear',      # Use linear if IDPP fails
        'options': {
            'fmax': 0.1,
            'steps': 100
        }
    }
]
```

### ⚙️ Global Settings

```python
pathway = build_reaction_pathway(
    steps=my_steps,
    base_indices=[0, 1, 2, 3, 4],
    global_settings={
        'default_method': 'idpp',     # Default interpolation method
        'default_frames': 15,         # Default frame count
        'default_fmax': 0.1          # Default IDPP tolerance
    }
)
```

---

## 📊 Performance Tips

### 🚀 Speed Optimization

1. **Use `'linear'` for simple transitions:**
   ```python
   'method': 'linear'  # ~1000x faster than IDPP
   ```

2. **Reduce frames for prototyping:**
   ```python
   'frames': 5  # Quick testing, increase later
   ```

3. **Skip frames in long trajectories:**
   ```python
   'skip': 5  # Use every 5th frame to reduce size
   ```

### 🎯 Quality Optimization

1. **Use `'idpp'` for chemical accuracy:**
   ```python
   'method': 'idpp'
   'options': {'fmax': 0.05}  # Higher precision
   ```

2. **Increase frames for smooth paths:**
   ```python
   'frames': 50  # More frames = smoother animation
   ```

3. **Use appropriate base_indices:**
   ```python
   base_indices=list(range(20))  # Include key atoms only
   ```

---

## 🚨 Common Issues and Solutions

### ❌ Issue: "Images have different numbers of atoms"
**Solution:** Use only atoms present in all structures for base_indices
```python
# Bad: using all atoms when structures have different sizes
base_indices = list(range(100))  

# Good: use common atoms only  
base_indices = list(range(50))   # Only first 50 atoms
```

### ❌ Issue: IDPP interpolation fails
**Solution:** Use fallback mechanism
```python
{
    'type': 'interpolate',
    'method': 'idpp',
    'fallback': 'linear',        # Auto-fallback to linear
    'options': {'fmax': 0.2}     # Relax convergence
}
```

### ❌ Issue: "n_frames must be at least 2"
**Solution:** Check frame count values
```python
'frames': 5     # ✅ Good: at least 2
'frames': 1     # ❌ Bad: too few frames
'frames': 0     # ❌ Bad: invalid
```

---

## 📁 Output and Visualization

### 💾 File Formats

The pathway is automatically saved in multiple formats:
```python
output_path='my_pathway.extxyz'  # Will create:
# - my_pathway.extxyz (extended XYZ with metadata)
# - my_pathway.traj (ASE trajectory format)
```

### 👁️ Visualization

```bash
# ASE GUI (recommended)
ase gui my_pathway.extxyz

# VMD
vmd my_pathway.extxyz

# Ovito
ovito my_pathway.extxyz

# PyMOL (convert first)
ase convert my_pathway.extxyz my_pathway.pdb
```

---

## 🎯 Best Practices Summary

1. **Choose the right step type:**
   - `interpolate`: For transitions between structures
   - `trajectory`: For existing MD/simulation data  
   - `frame`: For checkpoints and static states

2. **Method selection:**
   - `'idpp'`: Chemical accuracy, slower
   - `'linear'`: Speed, less accurate

3. **Frame management:**
   - More frames = smoother but larger files
   - Use `skip` parameter for large trajectories
   - Use `'auto'` frames for adaptive quality

4. **Error prevention:**
   - Always include fallback methods
   - Use appropriate base_indices
   - Test with small frame counts first

5. **Performance:**
   - Start with linear method for prototyping
   - Use IDPP for final high-quality results
   - Monitor progress with callback functions

Happy molecular pathway building! 🚀🧬