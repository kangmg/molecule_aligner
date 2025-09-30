# 🚀 Molecule Aligner Advanced Usage Guide

## 개요

Molecule Aligner v0.2.0은 단일 이미지, 궤적, IDPP 보간을 자유롭게 조합할 수 있는 통합 API를 제공합니다. 하나의 함수 호출로 복잡한 분자 동역학 경로를 생성할 수 있습니다.

## 🎯 핵심 철학

**"모든 것을 하나의 함수로"** - 단일 이미지, 궤적 파일, 보간 설정을 자유롭게 혼합하여 복잡한 반응 경로를 간단하게 구성

## 📚 통합 입력 형식

### 기본 구조
```python
reactions = [
    {
        # === 입력 소스 (택일) ===
        'traj_path': 'file.traj',           # 궤적 파일
        'traj': [atoms1, atoms2, ...],      # 직접 궤적
        'single_image': atoms,              # 단일 구조
        'single_image_path': 'struct.xyz',  # 단일 구조 파일
        
        # === 보간 설정 (single_image용) ===
        'interpolate_to': atoms,            # 보간 대상 구조
        'interpolate_to_path': 'target.xyz', # 보간 대상 파일
        'n_frames': 20,                     # 보간 프레임 수
        'interpolation_method': 'idpp',     # 'linear' or 'idpp'
        'interpolation_options': {          # 세부 옵션
            'fmax': 0.1,
            'steps': 100,
            'optimizer': 'LBFGS'
        },
        
        # === 공통 설정 ===
        'reverse': False,                   # 궤적 뒤집기
        'base_indices': [0, 1, 2, ...],    # 정렬 기준 원자
        
        # === 고급 설정 ===
        'skip_frames': 2,                   # 프레임 건너뛰기
        'frame_range': [10, 50],           # 특정 프레임 범위만 사용
        'weight': 1.0                       # 이 segment의 가중치 (미래 기능)
    }
]
```

## 🎮 사용 예시 모음

### 예시 1: 기본 IDPP 보간
```python
from molecule_aligner import align_and_merge_reactions_enhanced

# 두 구조 사이의 기본 IDPP 보간
reactions = [
    {
        'single_image_path': 'reactant.xyz',
        'interpolate_to_path': 'product.xyz',
        'n_frames': 15,
        'interpolation_method': 'idpp',
        'base_indices': [0, 1, 2, 3, 4, 5]
    }
]

trajectory = align_and_merge_reactions_enhanced(
    reactions=reactions,
    output_path='reaction_path.extxyz'
)
```

### 예시 2: 다중 보간 세그먼트 + 궤적 혼합
```python
# 복잡한 반응 경로: 보간 + 기존 궤적 + 보간
reactions = [
    # 세그먼트 1: 초기 접근 (IDPP)
    {
        'single_image': reactant_structure,
        'interpolate_to': intermediate1,
        'n_frames': 12,
        'interpolation_method': 'idpp',
        'interpolation_options': {
            'fmax': 0.05,        # 높은 정확도
            'optimizer': 'LBFGS',
            'steps': 150
        },
        'base_indices': [0, 1, 2, 3, 4]
    },
    
    # 세그먼트 2: 기존 MD 궤적 사용
    {
        'traj_path': 'md_simulation.traj',
        'reverse': False,
        'frame_range': [100, 200],  # 100-200번 프레임만 사용
        'skip_frames': 2,           # 매 2프레임마다 하나씩
        'base_indices': [0, 1, 2, 3, 4]
    },
    
    # 세그먼트 3: 최종 제품 형성 (선형 보간)
    {
        'single_image': intermediate2,
        'interpolate_to_path': 'final_product.xyz',
        'n_frames': 8,
        'interpolation_method': 'linear',  # 빠른 선형 보간
        'base_indices': [0, 1, 2, 3, 4]
    }
]

complete_trajectory = align_and_merge_reactions_enhanced(
    reactions=reactions,
    output_path='complex_reaction_pathway.extxyz',
    reference='first'
)
```

### 예시 3: 실제 촉매 반응 모델링
```python
# Suzuki 커플링 반응의 전체 사이클
suzuki_cycle = [
    # 1. 산화적 첨가 (Oxidative Addition)
    {
        'single_image_path': 'pd_catalyst.xyz',
        'interpolate_to_path': 'pd_aryl_halide.xyz', 
        'n_frames': 20,
        'interpolation_method': 'idpp',
        'interpolation_options': {'fmax': 0.02},  # 매우 정밀
        'base_indices': list(range(15))  # Pd 중심 + 리간드
    },
    
    # 2. 트랜스메탈레이션 전 과정 (MD 궤적)
    {
        'traj_path': 'transmetalation_md.traj',
        'reverse': False,
        'skip_frames': 5,  # 긴 MD에서 샘플링
        'base_indices': list(range(15))
    },
    
    # 3. 환원적 제거 (Reductive Elimination) - 고해상도
    {
        'single_image_path': 'pre_elimination.xyz',
        'interpolate_to_path': 'product_complex.xyz',
        'n_frames': 30,  # 많은 프레임으로 세밀하게
        'interpolation_method': 'idpp',
        'base_indices': list(range(15))
    },
    
    # 4. 제품 해리 - 빠른 과정
    {
        'single_image_path': 'product_complex.xyz', 
        'interpolate_to_path': 'regenerated_catalyst.xyz',
        'n_frames': 10,
        'interpolation_method': 'linear',  # 단순한 해리
        'base_indices': list(range(15))
    }
]

suzuki_trajectory = align_and_merge_reactions_enhanced(
    reactions=suzuki_cycle,
    output_path='suzuki_complete_cycle.extxyz'
)

print(f"Suzuki 사이클: {len(suzuki_trajectory)} 프레임 생성")
```

### 예시 4: 단백질 폴딩 경로
```python
# 단백질 펼침 → 중간체 → 최종 폴딩
protein_folding = [
    # 초기 펼침 상태 → 첫 번째 중간체 (느린 과정)
    {
        'single_image_path': 'unfolded_protein.pdb',
        'interpolate_to_path': 'intermediate1.pdb', 
        'n_frames': 50,  # 복잡한 구조 변화
        'interpolation_method': 'idpp',
        'interpolation_options': {
            'fmax': 0.15,     # 단백질은 유연하므로 관대하게
            'steps': 200,
            'optimizer': 'FIRE'  # 단백질에 적합한 optimizer
        },
        'base_indices': list(range(20, 35))  # 보존된 코어 영역만
    },
    
    # MD 시뮬레이션 결과 (중간체 동역학)
    {
        'traj_path': 'intermediate_dynamics.traj',
        'reverse': False,
        'frame_range': [0, 100],
        'skip_frames': 3,
        'base_indices': list(range(20, 35))
    },
    
    # 최종 폴딩 (빠른 붕괴)
    {
        'single_image_path': 'intermediate2.pdb',
        'interpolate_to_path': 'native_structure.pdb',
        'n_frames': 25,
        'interpolation_method': 'idpp', 
        'base_indices': list(range(20, 35))
    }
]
```

### 예시 5: 리간드 결합 과정의 상세 모델링
```python
# 약물-단백질 결합의 전 과정
drug_binding = [
    # 1. 원거리 접근 (긴 거리 → 단백질 표면)
    {
        'single_image_path': 'ligand_distant.xyz',
        'interpolate_to_path': 'ligand_surface.xyz',
        'n_frames': 15,
        'interpolation_method': 'linear',  # 단순한 확산
        'base_indices': list(range(50))    # 단백질 활성부위
    },
    
    # 2. 표면 탐색 및 인식 (MD 궤적)
    {
        'traj_path': 'surface_exploration.traj',
        'reverse': False,
        'skip_frames': 10,  # 긴 MD 시뮬레이션 샘플링
        'base_indices': list(range(50))
    },
    
    # 3. 포켓 진입 (유도 적합 모델)
    {
        'single_image_path': 'ligand_entrance.xyz',
        'interpolate_to_path': 'ligand_pocket.xyz', 
        'n_frames': 35,  # 복잡한 유도 적합
        'interpolation_method': 'idpp',
        'interpolation_options': {
            'fmax': 0.08,
            'steps': 120,
            'optimizer': 'LBFGS'
        },
        'base_indices': list(range(50))
    },
    
    # 4. 최적화 및 안정화
    {
        'traj_path': 'binding_optimization.traj',
        'reverse': False,
        'base_indices': list(range(50))
    }
]
```

### 예시 6: 역반응 포함 전체 사이클
```python
# 전방향 + 역방향 반응의 완전한 사이클
complete_cycle = [
    # Forward reaction: A → B
    {
        'single_image_path': 'reactant_A.xyz',
        'interpolate_to_path': 'product_B.xyz',
        'n_frames': 20,
        'interpolation_method': 'idpp',
        'base_indices': [0, 1, 2, 3, 4]
    },
    
    # 평형 상태에서의 동역학
    {
        'traj_path': 'equilibrium_dynamics.traj', 
        'reverse': False,
        'frame_range': [50, 150],
        'base_indices': [0, 1, 2, 3, 4]
    },
    
    # Reverse reaction: B → A (궤적 뒤집기 활용)
    {
        'single_image_path': 'product_B.xyz',
        'interpolate_to_path': 'reactant_A.xyz',
        'n_frames': 18,
        'interpolation_method': 'idpp',
        'reverse': True,  # 이 세그먼트를 뒤집어서 B→A로
        'base_indices': [0, 1, 2, 3, 4]
    }
]
```

## 🔧 고급 기능 활용

### 프레임 제어
```python
# 정밀한 프레임 제어 예시
reactions = [
    {
        'traj_path': 'long_simulation.traj',
        'frame_range': [100, 500],    # 100-500번 프레임만
        'skip_frames': 5,             # 매 5프레임마다 하나씩  
        'reverse': False,
        'base_indices': [0, 1, 2]
    },
    {
        'single_image_path': 'checkpoint.xyz',
        'interpolate_to_path': 'target.xyz',
        'n_frames': 25,               # 정확히 25개 중간 프레임
        'interpolation_method': 'idpp',
        'base_indices': [0, 1, 2] 
    }
]
```

### 적응형 보간 프레임 수
```python
def adaptive_frames(start_structure, end_structure, base_indices):
    """구조 차이에 따라 프레임 수를 자동 결정"""
    import numpy as np
    
    pos1 = start_structure.positions[base_indices]
    pos2 = end_structure.positions[base_indices]
    rmsd = np.sqrt(np.mean((pos1 - pos2)**2))
    
    if rmsd < 0.5:
        return 8   # 작은 변화
    elif rmsd < 2.0:
        return 15  # 중간 변화
    else:
        return 25  # 큰 변화

# 사용 예
reactions = [
    {
        'single_image': start_struct,
        'interpolate_to': end_struct, 
        'n_frames': adaptive_frames(start_struct, end_struct, [0,1,2,3]),
        'interpolation_method': 'idpp',
        'base_indices': [0, 1, 2, 3]
    }
]
```

### 조건부 IDPP 사용
```python
# IDPP 실패시 자동 linear 전환
reactions = [
    {
        'single_image_path': 'complex_start.xyz',
        'interpolate_to_path': 'complex_end.xyz',
        'n_frames': 20,
        'interpolation_method': 'idpp',
        'interpolation_options': {
            'fmax': 0.1,
            'steps': 100,
            'fallback_to_linear': True  # IDPP 실패시 linear로 자동 전환
        },
        'base_indices': list(range(20))
    }
]
```

## 📊 품질 모니터링

### 궤적 품질 실시간 분석
```python
from molecule_aligner.interpolate import interpolation_quality_check

# 궤적 생성 후 품질 체크
trajectory = align_and_merge_reactions_enhanced(reactions)

# 각 세그먼트별 품질 분석
segment_starts = [0, 20, 45, 70]  # 각 세그먼트 시작점
for i in range(len(segment_starts) - 1):
    start_idx = segment_starts[i] 
    end_idx = segment_starts[i + 1]
    segment = trajectory[start_idx:end_idx]
    
    if len(segment) >= 2:
        quality = interpolation_quality_check(
            segment, 
            segment[0], segment[-1],
            f"Segment_{i+1}"
        )
        print(f"세그먼트 {i+1} 품질:")
        print(f"  - 경로 길이: {quality['total_path_length']:.2f} Å")
        print(f"  - 부드러움: {quality['smoothness_score']:.3f}")
        print(f"  - 평균 스텝: {quality['avg_step_size']:.3f} Å")
```

## 🎯 성능 최적화 팁

### 1. 보간 방법 선택 가이드
```python
# 빠른 프로토타이핑
'interpolation_method': 'linear'

# 화학적 정확성이 중요한 경우  
'interpolation_method': 'idpp'
'interpolation_options': {'fmax': 0.05, 'steps': 150}

# 큰 시스템이나 성능이 중요한 경우
'interpolation_method': 'idpp' 
'interpolation_options': {'fmax': 0.2, 'steps': 50}
```

### 2. 메모리 효율적 처리
```python
# 큰 궤적의 경우 청크 단위 처리
def process_large_trajectory(file_path, chunk_size=100):
    reactions = []
    
    # 궤적을 청크로 나누어 처리
    total_frames = len(read(file_path, ':'))
    
    for start in range(0, total_frames, chunk_size):
        end = min(start + chunk_size, total_frames)
        
        reactions.append({
            'traj_path': file_path,
            'frame_range': [start, end],
            'skip_frames': 2,  # 메모리 절약
            'base_indices': [0, 1, 2, 3, 4]
        })
    
    return reactions
```

### 3. 병렬 처리 준비
```python
# 독립적인 세그먼트들을 병렬로 처리 준비
parallel_segments = [
    # 각 세그먼트가 독립적으로 처리 가능
    [segment1_reactions],
    [segment2_reactions], 
    [segment3_reactions]
]

# 나중에 합치기
# final_trajectory = merge_parallel_results(parallel_segments)
```

## 🚀 실전 워크플로우 예시

### 전체 약물 설계 파이프라인
```python
def drug_design_workflow(target_protein, candidate_ligands):
    """완전한 약물-표적 상호작용 분석"""
    
    all_trajectories = []
    
    for i, ligand in enumerate(candidate_ligands):
        binding_pathway = [
            # 1. 접근 단계
            {
                'single_image_path': f'ligand_{i}_distant.xyz',
                'interpolate_to_path': f'ligand_{i}_surface.xyz',
                'n_frames': 15,
                'interpolation_method': 'linear',
                'base_indices': list(range(50))
            },
            
            # 2. 결합 단계  
            {
                'single_image_path': f'ligand_{i}_surface.xyz',
                'interpolate_to_path': f'ligand_{i}_bound.xyz',
                'n_frames': 30,
                'interpolation_method': 'idpp',
                'interpolation_options': {'fmax': 0.1},
                'base_indices': list(range(50))
            },
            
            # 3. 최적화 단계 (MD 결과)
            {
                'traj_path': f'optimization_{i}.traj',
                'skip_frames': 5,
                'base_indices': list(range(50))
            }
        ]
        
        # 각 리간드별 궤적 생성
        trajectory = align_and_merge_reactions_enhanced(
            reactions=binding_pathway,
            output_path=f'ligand_{i}_binding_pathway.extxyz'
        )
        
        all_trajectories.append(trajectory)
    
    return all_trajectories
```

## 🎓 모범 사례

### 1. 구조화된 프로젝트 조직
```python
# 프로젝트 구조
project/
├── structures/           # 입력 구조들
├── trajectories/        # MD 궤적들  
├── analysis/           # 분석 스크립트들
└── results/            # 결과 궤적들

# 설정 파일 사용
import yaml

with open('reaction_config.yaml') as f:
    config = yaml.load(f)

reactions = config['reaction_steps']
trajectory = align_and_merge_reactions_enhanced(reactions)
```

### 2. 오류 처리 및 로깅
```python
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def robust_trajectory_generation(reactions):
    """견고한 궤적 생성 with 오류 처리"""
    
    for i, reaction in enumerate(reactions):
        try:
            # 각 반응 단계 검증
            validate_reaction_step(reaction)
            logger.info(f"단계 {i+1} 검증 완료")
            
        except Exception as e:
            logger.error(f"단계 {i+1} 오류: {e}")
            # 대체 방법 시도
            reaction['interpolation_method'] = 'linear'
            logger.info(f"단계 {i+1} linear 방법으로 대체")
    
    return align_and_merge_reactions_enhanced(reactions)
```

### 3. 재현 가능한 연구
```python
# 시드 설정으로 재현 가능한 결과
import numpy as np
np.random.seed(42)

# 모든 매개변수 문서화
trajectory_metadata = {
    'creation_date': '2024-01-01',
    'method': 'IDPP_interpolation',
    'parameters': {
        'fmax': 0.1,
        'optimizer': 'LBFGS',
        'base_indices': [0, 1, 2, 3, 4]
    },
    'input_files': ['reactant.xyz', 'product.xyz'],
    'total_frames': len(trajectory)
}

# 메타데이터와 함께 저장
with open('trajectory_metadata.json', 'w') as f:
    json.dump(trajectory_metadata, f, indent=2)
```

이 가이드를 통해 Molecule Aligner의 모든 기능을 최대한 활용하여 복잡한 분자 동역학 연구를 효율적으로 수행할 수 있습니다! 🚀