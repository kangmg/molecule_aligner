# 🎯 통합 API 제안: Universal Reaction Builder

## 목표

사용자가 요청한 "single image, traj(reverse:bool) 등을 섞어서 interpolate의 경우 각각 수를 지정해서 하나의 함수로 실행" 할 수 있는 통합 API 설계

## 🚀 새로운 통합 함수: `build_reaction_pathway()`

```python
from molecule_aligner import build_reaction_pathway

def build_reaction_pathway(
    steps: List[Dict],
    base_indices: List[int],
    output_path: Optional[str] = None,
    reference: str = 'first',
    global_settings: Optional[Dict] = None
) -> List[Atoms]:
    """
    통합 반응 경로 생성기 - 모든 유형의 입력을 하나의 함수로 처리
    
    Args:
        steps: 반응 단계들의 리스트
        base_indices: 전역 정렬 기준 원자 인덱스
        output_path: 결과 저장 경로
        reference: 기준 프레임 ('first', 'reactant')
        global_settings: 전역 설정 (기본값들)
    
    Returns:
        완성된 반응 경로 (List[Atoms])
    """
```

## 📋 통합 입력 형식

### Step Type 1: Single Image with Interpolation
```python
{
    'type': 'interpolate',
    'from': 'start.xyz',              # 시작 구조 (파일 경로 또는 Atoms)
    'to': 'end.xyz',                  # 끝 구조  
    'frames': 20,                     # 보간 프레임 수
    'method': 'idpp',                 # 'linear' or 'idpp'
    'options': {                      # 메소드별 세부 옵션
        'fmax': 0.1,
        'optimizer': 'LBFGS',
        'steps': 100
    },
    'reverse': False                   # 결과 뒤집기
}
```

### Step Type 2: Trajectory Segment  
```python
{
    'type': 'trajectory',
    'source': 'simulation.traj',       # 궤적 파일 또는 List[Atoms]
    'frames': 'all',                   # 'all', [start, end], 또는 특정 수
    'skip': 1,                         # 프레임 건너뛰기 (1=모두, 2=하나걸러)  
    'reverse': True                    # 궤적 뒤집기
}
```

### Step Type 3: Single Frame Insert
```python
{
    'type': 'frame',
    'source': 'intermediate.xyz',      # 단일 구조
    'repeat': 3                        # 동일 프레임 반복 횟수 (선택)
}
```

## 🎮 실제 사용 예시

### 예시 1: 기본 IDPP 보간 체인
```python
# A → B → C 연쇄 보간
pathway = build_reaction_pathway(
    steps=[
        {
            'type': 'interpolate',
            'from': 'A.xyz',
            'to': 'B.xyz', 
            'frames': 15,
            'method': 'idpp'
        },
        {
            'type': 'interpolate', 
            'from': 'B.xyz',
            'to': 'C.xyz',
            'frames': 12,
            'method': 'idpp'
        }
    ],
    base_indices=[0, 1, 2, 3, 4],
    output_path='A_to_C_pathway.extxyz'
)
```

### 예시 2: 복잡한 혼합 워크플로우
```python
# 실제 복잡한 반응: 보간 + MD궤적 + 역방향 + 보간
complex_reaction = build_reaction_pathway(
    steps=[
        # 1. 초기 접근 (IDPP 보간)
        {
            'type': 'interpolate',
            'from': reactant_atoms,
            'to': 'intermediate1.xyz',
            'frames': 25,
            'method': 'idpp',
            'options': {'fmax': 0.05, 'optimizer': 'LBFGS'}
        },
        
        # 2. MD 시뮬레이션 결과 사용
        {
            'type': 'trajectory',
            'source': 'md_run1.traj',
            'frames': [50, 150],  # 50-150번 프레임만
            'skip': 3,            # 매 3번째만 선택
            'reverse': False
        },
        
        # 3. 중간 안정화 지점
        {
            'type': 'frame',
            'source': 'stable_intermediate.xyz',
            'repeat': 5  # 5프레임 동안 정적 상태
        },
        
        # 4. 역반응 경로 (기존 궤적 뒤집기)
        {
            'type': 'trajectory',
            'source': forward_trajectory,  # List[Atoms]
            'frames': 30,     # 처음 30개만
            'reverse': True   # 뒤집어서 역반응으로
        },
        
        # 5. 최종 제품 형성 (빠른 선형)
        {
            'type': 'interpolate',
            'from': 'pre_product.xyz',
            'to': 'final_product.xyz',
            'frames': 8,
            'method': 'linear'  # 빠른 마무리
        }
    ],
    base_indices=list(range(20)),  # 처음 20개 원자로 정렬
    output_path='complete_reaction_cycle.extxyz'
)

print(f"완성된 반응 경로: {len(complex_reaction)} 프레임")
```

### 예시 3: 사용자 요청과 정확히 일치하는 예시
```python
# "single image, traj(reverse:bool) 등을 섞어서 interpolate의 경우 각각 수를 지정"
user_workflow = build_reaction_pathway(
    steps=[
        # Single image → interpolate (사용자 지정 프레임 수)
        {
            'type': 'interpolate',
            'from': 'reactant.xyz',      # single image
            'to': 'ts.xyz',              # single image  
            'frames': 18,                # 사용자 지정 수
            'method': 'idpp'
        },
        
        # Trajectory with reverse
        {
            'type': 'trajectory',
            'source': 'ts_to_product.traj',  # traj
            'reverse': True,                  # reverse: bool
            'frames': 25                      # 사용자 지정 수
        },
        
        # Another interpolate with different frame count
        {
            'type': 'interpolate', 
            'from': 'product1.xyz',     # single image
            'to': 'product2.xyz',       # single image
            'frames': 12,               # 다른 프레임 수
            'method': 'linear'
        },
        
        # Mixed: single frame + trajectory + interpolate
        {
            'type': 'frame',
            'source': 'catalyst.xyz'    # single frame insert
        },
        {
            'type': 'trajectory',
            'source': existing_traj,    # traj (List[Atoms])
            'reverse': False,           # reverse: bool
            'frames': 20
        },
        {
            'type': 'interpolate',
            'from': 'intermediate.xyz',
            'to': 'final.xyz', 
            'frames': 15,               # 또 다른 사용자 지정 수
            'method': 'idpp'
        }
    ],
    base_indices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20],  # 0-20번 원자
    output_path='user_custom_pathway.extxyz'
)
```

## 🔧 고급 기능

### 조건부 처리
```python
# IDPP 실패시 자동 linear 전환
steps = [
    {
        'type': 'interpolate',
        'from': 'complex_start.xyz',
        'to': 'complex_end.xyz',
        'frames': 20,
        'method': 'idpp',
        'fallback': 'linear',        # 실패시 대체 방법
        'options': {'fmax': 0.1, 'steps': 100}
    }
]
```

### 적응형 프레임 수
```python
# 구조 차이에 따른 자동 프레임 수 결정
steps = [
    {
        'type': 'interpolate',
        'from': 'A.xyz',
        'to': 'B.xyz', 
        'frames': 'auto',           # 자동 결정
        'frames_per_angstrom': 5,   # Å당 5프레임
        'method': 'idpp'
    }
]
```

### 품질 기반 적응
```python
# 품질 모니터링 및 자동 조정
steps = [
    {
        'type': 'interpolate',
        'from': 'start.xyz',
        'to': 'end.xyz',
        'frames': 15,
        'method': 'idpp',
        'quality_check': True,       # 품질 검사 활성화
        'min_quality_score': 0.8,    # 최소 품질 요구사항
        'auto_refine': True          # 품질 부족시 자동 개선
    }
]
```

## 📊 실시간 진행 모니터링

```python
# 진행상황 콜백 함수
def progress_callback(step_index, step_type, progress_percent):
    print(f"Step {step_index} ({step_type}): {progress_percent}% 완료")

pathway = build_reaction_pathway(
    steps=my_steps,
    base_indices=my_indices,
    progress_callback=progress_callback,  # 진행상황 모니터링
    output_path='monitored_pathway.extxyz'
)
```

## 🎯 구현 계획

### Phase 1: 핵심 통합 함수
```python
def build_reaction_pathway(steps, base_indices, **kwargs):
    """핵심 통합 함수 구현"""
    
    all_reactions = []
    
    for step in steps:
        if step['type'] == 'interpolate':
            reaction = _handle_interpolate_step(step, base_indices)
        elif step['type'] == 'trajectory': 
            reaction = _handle_trajectory_step(step, base_indices)
        elif step['type'] == 'frame':
            reaction = _handle_frame_step(step, base_indices)
        
        all_reactions.append(reaction)
    
    return align_and_merge_reactions_enhanced(
        reactions=all_reactions,
        **kwargs
    )
```

### Phase 2: 고급 기능 추가
- 적응형 프레임 수 결정
- 품질 기반 자동 조정  
- 진행상황 모니터링
- 병렬 처리 지원

### Phase 3: 사용자 편의 기능
- 설정 파일 지원 (YAML/JSON)
- 시각화 도구 통합
- 템플릿 라이브러리

## 💡 사용자 요청 완벽 구현

```python
# 사용자 요청: "single image, traj(reverse:bool) 등을 섞어서 
# interpolate의 경우 각각 수를 지정해서 하나의 함수로 실행"

from molecule_aligner import build_reaction_pathway

# ✅ 정확히 요청된 기능
mixed_pathway = build_reaction_pathway(
    steps=[
        # single image → interpolate (프레임 수 개별 지정)
        {'type': 'interpolate', 'from': 'img1.xyz', 'to': 'img2.xyz', 'frames': 15},
        
        # traj with reverse=True
        {'type': 'trajectory', 'source': 'traj1.traj', 'reverse': True, 'frames': 20},
        
        # another interpolate (다른 프레임 수)
        {'type': 'interpolate', 'from': 'img3.xyz', 'to': 'img4.xyz', 'frames': 8},
        
        # traj with reverse=False  
        {'type': 'trajectory', 'source': 'traj2.traj', 'reverse': False, 'frames': 12},
        
        # final interpolate (또 다른 프레임 수)
        {'type': 'interpolate', 'from': 'img5.xyz', 'to': 'img6.xyz', 'frames': 25}
    ],
    base_indices=list(range(71)),  # 0-70 원자로 정렬
    output_path='user_requested_pathway.extxyz'
)

# ✅ 하나의 함수로 모든 것을 처리!
print(f"완성: {len(mixed_pathway)} 프레임")
```

이 통합 API는 사용자가 요청한 모든 기능을 하나의 직관적인 함수로 제공하며, 확장성과 유연성을 동시에 보장합니다! 🚀