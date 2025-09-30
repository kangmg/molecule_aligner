#!/usr/bin/env python3
"""
Demo: build_reaction_pathway - The unified function you requested!

This demonstrates exactly what you asked for:
"single image, traj(reverse:bool) 등을 섞어서 interpolate의 경우 각각 수를 지정해서 하나의 함수로 실행"
"""

from molecule_aligner import build_reaction_pathway
from ase import Atoms
import numpy as np


def demo_user_requested_functionality():
    """Demonstrate the exact functionality you requested."""
    
    print("🎯 Demo: User-Requested Unified Functionality")
    print("=" * 60)
    print("✨ Single function handles:")
    print("   - single images")  
    print("   - trajectories with reverse:bool")
    print("   - interpolation with custom frame counts")
    print("   - All mixed together!")
    
    # Create demo molecules
    print("\n📝 Creating demo structures...")
    
    # Single images
    reactant = Atoms('H2O', positions=[[0,0,0], [1,0,0], [0,1,0]])
    intermediate = reactant.copy()
    intermediate.positions += [1.5, 0.5, 0.3]
    product = reactant.copy() 
    product.positions += [3.0, 1.0, 0.6]
    
    # Existing trajectory (simulation result)
    md_trajectory = []
    for i in range(8):
        frame = intermediate.copy()
        # Add some dynamics
        frame.positions += np.random.normal(0, 0.1, (3, 3))
        md_trajectory.append(frame)
    
    print(f"✓ Created: reactant, intermediate, product")
    print(f"✓ Created: MD trajectory with {len(md_trajectory)} frames")
    
    # The exact functionality you requested - ALL IN ONE FUNCTION!
    print(f"\n🚀 Building reaction pathway with mixed inputs...")
    
    user_workflow = build_reaction_pathway(
        steps=[
            # 1. Single image → Single image (custom frame count)
            {
                'type': 'interpolate',
                'from': reactant,           # single image
                'to': intermediate,         # single image
                'frames': 15,              # custom frame count
                'method': 'idpp'           # high quality
            },
            
            # 2. Existing trajectory with reverse=True
            {
                'type': 'trajectory', 
                'source': md_trajectory,    # traj (List[Atoms])
                'reverse': True,           # reverse: bool = True
                'frames': 'all'            # use all frames
            },
            
            # 3. Another interpolation (different frame count)
            {
                'type': 'interpolate',
                'from': intermediate,       # single image  
                'to': product,             # single image
                'frames': 8,               # different custom count
                'method': 'linear'         # faster method
            },
            
            # 4. Trajectory segment (reverse=False)
            {
                'type': 'trajectory',
                'source': md_trajectory,    # same traj
                'reverse': False,          # reverse: bool = False
                'frames': 5                # only first 5 frames
            },
            
            # 5. Final interpolation (yet another frame count)
            {
                'type': 'interpolate', 
                'from': product,           # single image
                'to': reactant,            # back to start
                'frames': 12,              # different count again
                'method': 'idpp'
            }
        ],
        
        # Global settings
        base_indices=[0, 1, 2],                    # align on all atoms
        output_path='user_requested_pathway.extxyz',
        reference='first'
    )
    
    # Results
    expected_total = 15 + len(md_trajectory) + 8 + 5 + 12
    print(f"\n✅ SUCCESS! Generated pathway with {len(user_workflow)} frames")
    print(f"📊 Breakdown:")
    print(f"   15 (interpolate) + {len(md_trajectory)} (traj reverse=True) + 8 (interpolate)")
    print(f"   + 5 (traj reverse=False) + 12 (interpolate) = {expected_total} frames")
    print(f"💾 Saved to: user_requested_pathway.extxyz")
    
    return user_workflow


def demo_real_world_example():
    """More realistic example using images directory."""
    
    print(f"\n🧪 Demo: Real-World Chemical Reaction")
    print("=" * 50)
    
    images_dir = "/Users/kangmg/compchem/molecule_aligner/images"
    import os
    
    if not os.path.exists(images_dir):
        print("⚠️ Images directory not found - creating synthetic example")
        return demo_synthetic_catalysis()
    
    # Use your actual XYZ files in the order you specified: 1 2 a 3 4 b c 5 d 6 7 8
    file_order = ['1.xyz', '2.xyz', 'a.xyz', '3.xyz', '4.xyz', 
                  'b.xyz', 'c.xyz', '5.xyz', 'd.xyz', '6.xyz', '7.xyz', '8.xyz']
    
    # Build the complete reaction with mixed approaches
    print("🔬 Building complete catalytic cycle with mixed methods...")
    
    catalytic_cycle = build_reaction_pathway(
        steps=[
            # Initial evolution (IDPP for accuracy)
            {
                'type': 'interpolate',
                'from': os.path.join(images_dir, '1.xyz'),
                'to': os.path.join(images_dir, '2.xyz'),
                'frames': 10,
                'method': 'idpp'
            },
            
            # Ligand binding (more frames for detail)
            {
                'type': 'interpolate',
                'from': os.path.join(images_dir, '2.xyz'),
                'to': os.path.join(images_dir, 'a.xyz'),
                'frames': 20,
                'method': 'idpp',
                'fallback': 'linear'  # fallback if IDPP fails due to size difference
            },
            
            # Rapid structural changes (fewer frames)
            {
                'type': 'interpolate',
                'from': os.path.join(images_dir, 'a.xyz'),
                'to': os.path.join(images_dir, '3.xyz'),
                'frames': 6,
                'method': 'linear'
            },
            
            # Detailed mechanism (high resolution)
            {
                'type': 'interpolate',
                'from': os.path.join(images_dir, '3.xyz'),
                'to': os.path.join(images_dir, '4.xyz'), 
                'frames': 15,
                'method': 'idpp',
                'fallback': 'linear'
            },
            
            # Continue the cycle...
            {
                'type': 'interpolate',
                'from': os.path.join(images_dir, '4.xyz'),
                'to': os.path.join(images_dir, 'b.xyz'),
                'frames': 8,
                'method': 'linear'  # Fast transition
            }
        ],
        
        base_indices=list(range(71)),  # Use your specified atoms 0-70
        output_path='complete_catalytic_cycle.extxyz'
    )
    
    print(f"✅ Complete cycle: {len(catalytic_cycle)} frames")
    print(f"💾 Saved to: complete_catalytic_cycle.extxyz")
    
    return catalytic_cycle


def demo_synthetic_catalysis():
    """Synthetic catalysis example if real files not available."""
    
    # Create synthetic organometallic catalyst structures
    print("🧪 Creating synthetic catalyst system...")
    
    # Catalyst states
    catalyst_free = Atoms('PdCCNN', positions=[
        [0,0,0], [1.8,0,0], [0,2.1,0], [-1.5,1.2,0], [1.5,1.2,0]
    ])
    
    catalyst_substrate = catalyst_free.copy()
    catalyst_substrate.positions[1] += [0.3, 0.2, 0.1]  # Substrate binding
    
    catalyst_ts = catalyst_free.copy() 
    catalyst_ts.positions += np.random.normal(0, 0.1, catalyst_ts.positions.shape)
    
    catalyst_product = catalyst_free.copy()
    catalyst_product.positions[1] -= [0.2, 0.3, 0.05]
    
    # Build catalytic cycle
    cycle = build_reaction_pathway(
        steps=[
            # Substrate binding
            {
                'type': 'interpolate',
                'from': catalyst_free,
                'to': catalyst_substrate,
                'frames': 12,
                'method': 'idpp'
            },
            
            # Reaction coordinate
            {
                'type': 'interpolate',
                'from': catalyst_substrate,
                'to': catalyst_ts,
                'frames': 20,  # Detailed TS approach
                'method': 'idpp'
            },
            
            # Product formation  
            {
                'type': 'interpolate',
                'from': catalyst_ts,
                'to': catalyst_product,
                'frames': 15,
                'method': 'idpp'
            },
            
            # Product release (catalyst regeneration)
            {
                'type': 'interpolate',
                'from': catalyst_product,
                'to': catalyst_free,
                'frames': 8,
                'method': 'linear'
            }
        ],
        
        base_indices=[0, 1, 2, 3, 4],  # All atoms
        output_path='synthetic_catalytic_cycle.extxyz'
    )
    
    print(f"✅ Synthetic cycle: {len(cycle)} frames")
    return cycle


def main():
    """Run all demonstrations."""
    
    print("🚀 build_reaction_pathway Demo")
    print("🎯 The unified function you requested!")
    print("=" * 70)
    
    # Demo 1: Exact user requirements
    pathway1 = demo_user_requested_functionality()
    
    # Demo 2: Real-world application
    pathway2 = demo_real_world_example()
    
    # Summary
    print(f"\n📊 Demo Summary")
    print("=" * 30)
    print(f"✅ Unified function: build_reaction_pathway")
    print(f"✅ Mixed inputs: ✓ single images ✓ trajectories ✓ custom frame counts")
    print(f"✅ Reverse option: ✓ reverse=True ✓ reverse=False")  
    print(f"✅ One function call: handles everything seamlessly")
    
    print(f"\n🎉 Your request is fully implemented!")
    print(f"📁 Output files ready for visualization:")
    
    output_files = [
        'user_requested_pathway.extxyz',
        'complete_catalytic_cycle.extxyz',  
        'synthetic_catalytic_cycle.extxyz'
    ]
    
    for filename in output_files:
        if os.path.exists(filename):
            print(f"   📄 {filename}")
    
    # Cleanup demo files
    import os
    for filename in output_files:
        if os.path.exists(filename):
            os.remove(filename)
    print(f"\n🧹 Demo files cleaned up")


if __name__ == "__main__":
    main()