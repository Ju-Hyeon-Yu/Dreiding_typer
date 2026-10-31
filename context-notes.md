# 작업 메모 (결정과 근거)

## 2026-10-26 출발점
- LUNAR DREIDING 타이핑은 말단 원자 혼성을 결합각으로 추정해 술포네이트 O가
  `O_31`(united-atom OH)로 들어가고, I/Br/P/Ca 분기가 코드 오타로 실패함 →
  자체 타이퍼로 결정.
- 사용자 결정: 출력은 `.lt`, RESP mol2는 SMILES로 결합 차수를 입혀 원자
  순서·전하 유지, 수소결합 항 포함.

## 기준 force field 파일
- moltemplate 배포본 `force_fields/dreiding.lt` (Matthew Bone, 2020).
  라벨 체계는 m-bone/moltemplate-DREIDING의 `DREIDING_Label_Manual.pdf`.
- 결합 차수는 라벨 쌍으로 결정됨 (생성 스크립트 `Dreiding_bond.py`의 if 순서):
  b1-b1 → 1, b2-b2 → 2, b2와 `_2` → 1, `_1`-`_1` → 3, `_2`-`_2` → 2,
  `_R`-`_R` → 1.5, 나머지 → 1.
- d 플래그 = 해당 원자의 이웃 수 − 1 (이면각 장벽 나눗셈에 쓰임).
- 허용 라벨은 실행 시 dreiding.lt의 `Data Masses`에서 읽어 버전 차이를 막음.

## 타입 규칙 결정
- 방향족(RDKit) → `_R`.
- 공명 등가 말단 O 묶음(카복실레이트, 니트로)의 중심 C/N과 O → `_R`
  (C_R–O_R 결합 차수 1.5로 표현 가능). 중심이 S/P면 DREIDING에 다중결합
  S/P가 없으므로 O는 `O_2`, 결합은 1로 둠(경고).
- sp2/방향족 원자에 붙은 3배위 비방향족 N(아마이드, 아닐린) → `N_R`.
  아마이드 C–N은 결합 차수 1, 이면각은 규칙 (e) V=5.
  (라벨 매뉴얼 예시는 아마이드 N을 `N_2_b1_d2_hd` + 카보닐 `C_2`로 적어
  C–N이 결합 차수 2가 됨 — 따르지 않음.)
- 공액 O(에스터, 페놀, 산 OH)는 `O_3` 유지. 원 논문 이면각 규칙 (i)
  (O족 sp3 + sp2/공명 원자, V=2, n=2)가 바로 이 경우를 위해 있고,
  dreiding.lt의 Case I가 이를 구현함.
- 다중결합 없는 2배위 N(이미드 음이온 등) → `N_R_d1` (경고 없이).
- S는 `S_3`만 존재. 3배위 이상 S(술포네이트·술폰)는 각도 92.1°,
  이면각 개수 가정이 실제와 다름 → 경고.
- Cl 질량이 dreiding.lt에 34.45로 들어가 있음(실제 35.45) → Cl 있으면 경고.

## 수소결합
- N/O/F에 붙은 H → `H_HB`, 그 무거운 원자 → `_hd`.
- 받개 `_ha`: 모든 O(주개 제외), F, 비공유전자쌍이 남는 N
  (방향족 2배위, sp3 3배위 이하 중성/음이온, 이민 N, 니트릴 N, 2배위 N_R).
- `_hd`/`_ha` 변형 라벨이 없는 조합은 경고 후 플래그 없이 둠.

## 테스트에서 드러난 것 (2026-10-26)
- 처음에는 결합 차수를 생성 스크립트 if 순서를 복제한 함수로 판정했는데,
  아닐라이드에서 `C_2`–`N_R_b1_d2_hd` 쌍이 결합 차수는 맞지만 dreiding.lt 에
  이면각 규칙이 없어 moltemplate 가 이면각 4개를 **경고 없이 뺐음**.
  → dreiding.lt 의 By Type 규칙(나중 규칙 우선)을 직접 읽어 결합 차수(K/350)와
  이면각 규칙 존재를 판정하도록 변경(DreidingFF). 결과: `C_2_b1`–`N_R_b1_d2_hd`.
- 이면각 계수를 Mayo 1990 규칙 (a)~(j)의 독립 구현(lammps_interface 이식)과
  대조한 결과, 차이는 S_3(3·4배위)와 4배위 N_3 의 이웃 수 가정뿐이었음 → 경고로 표시.
- SMILES만 줄 때 RDKit 이 이온쌍을 한 번에 임베딩하면 두 이온이 0.6 Å 까지
  겹쳐 LAMMPS LJ 에너지가 1e9 가 나옴 → 조각마다 따로 임베딩하고 x 축으로 4 Å 띄움.
- 진공 직접합: in.init/in.settings 의 lj/cut/coul/long → lj/cut/coul/cut,
  kspace_style 삭제, boundary f f f 로 LAMMPS 2024-08-29 에서 읽힘 확인.

## 2단계 결정 (2026-10-26)
- 사용자 요청: moltemplate.sh 까지 연결, 라벨은 dreiding.lt 것을 최대한 살리되 정확도 우선.
- 아마이드는 현재 방식 유지(C–N 단일, 장벽 5) — 매뉴얼 방식(C–N 이중, 장벽 45)보다
  실제(1.33 Å, 15–20 kcal/mol)에 가까움. 사용자가 "추천대로" 선택.
- 매뉴얼 예시와 라벨 이름이 다른 경우(사이클로펜타디에논, 이미다졸 N–H)는 이름 맞추기
  대신 moltemplate 결과 계수가 같음을 테스트로 보장. 매뉴얼은 예시만 있고 일반 규칙이
  아니라서 이름을 강제로 맞출 기준이 없음.
- cleanup_moltemplate.sh 로 타입 수가 107→8 (atom), 17576→5 (dihedral) 로 줄어듦.
  단, 주개가 없으면 hbond pair_coeff 가 모두 사라져 LAMMPS 가
  "hbond/dreiding/lj is not used" 로 멈춤 → pair_style 에서 제거해야 함.
  hybrid/overlay 에 sub-style 하나만 남겨도 LAMMPS 가 읽음(확인).

## 범위 확정 (2026-10-26)
- 사용자: 출력은 단일 분자 `.lt`(import "dreiding.lt" + `<이름> inherits DREIDING`)까지만.
  system.lt 와 moltemplate.sh 실행은 사용자가 직접 함 → `--lammps` 파이프라인 되돌림.
- moltemplate.sh 는 테스트에서 .lt 검증 용도로만 씀.

## 3단계 — 오프라인 HTML (2026-10-28)
- 요구: SMILES 필수 + (mol2/pdb/xyz 선택) → 타입 할당 .lt 다운로드·복사, 3D 뷰어에 타입 라벨,
  타입을 한눈에 보는 요약. 팀 배포용이라 회사 Edge·인터넷 연결 없이 단일 파일로 동작해야 함.
- 엔진: RDKit.js(MinimalLib 2026.09.1, Python RDKit 과 같은 버전) — 방향족성·AddHs 원자 순서가
  Python 과 같음. wasm 7.4 MB → gzip+base64 3.2 MB 로 내장, DecompressionStream 으로 풀기.
  Pyodide 는 내장 시 수십 MB 이고 이 환경에서 rdkit 휠을 받을 수 없어 제외.
- RDKit.js 에는 3D 생성·결합 추정이 없음 → SMILES 필수(사용자 결정), SMILES 만일 때는
  2D 좌표에서 DREIDING 결합·각 항 + 반발로 JS 최소화(사용자 결정 (a)).
- 전하: mol2 는 파일 값, pdb/xyz/SMILES 만은 0 + 경고 (사용자 결정).
- 규칙표는 빌드 때 Python DreidingFF 로 미리 계산해 JSON 으로 내장 → 판정이 Python 과 같음.
- JS 타이퍼는 Python 의 반복 순서(후보 순서, 결합 순서, BFS, 동률 처리)까지 그대로 옮겨
  결과가 같도록 함. mol2 입력일 때는 Python 처럼 원자를 파일 순서로 재번호한 뒤 타이핑.

### 3단계 구현 중 드러난 것
- 이 RDKit.js 빌드는 `wasmBinary` 옵션을 읽지 않아 file:// 에서 wasm 을 fetch 하다 실패함
  → `instantiateWasm` 콜백으로 압축 해제한 바이트를 직접 넘김.
- `.tip`, `.toolbar` 등 display 를 지정한 요소에서 `hidden` 속성이 무시됨 → `[hidden]{display:none!important}`.
- 3D 생성: 4배위와 3배위 이상 S·P 중심은 사면체(109.47°)를 목표로 둠(DREIDING S_3 92.1° 를 쓰면
  술포네이트가 찌그러짐). 3·4원 고리는 각이 맞지 않는 것이 정상이라 테스트에서 제외.
  SMILES 만으로 만든 이온쌍 .lt 를 moltemplate→LAMMPS 진공 최소화하면 PE 75.8 → 60.5 kcal/mol 로
  시작 구조가 무난함을 확인.
- 뷰어 색 기본값은 '타입' — 타입 카드 색과 1:1 로 맞아 한눈에 확인하기 쉬움(사용자 요청).
- 결과: 단일 HTML 3.57 MB, file:// 로딩 약 0.5 초, 외부 요청 0.

## 4단계 — 타입 정확도 검토 (2026-10-31)
- 기준: Caltech MSC(Goddard 그룹) dreid-typer v2.2 공식 규칙과 테스트 케이스(논문 CSD 분자·아미노산·핵산
  40개, 원자 730개), SEAMM dreiding.frc(독립 전사본)로 파라미터 대조. 원 논문 PDF 는 접근 불가.
- 우리 타이퍼와 다른 원자 65개(+타이핑 실패 4개 분자):
  카복실레이트·니트로 말단 O (우리 O_R, dreid-typer O_2), 아마이드 C (C_2 vs C_R),
  구아니디늄 (C_2/N_2 vs C_R/N_R), 술폰아마이드 N (N_R vs N_3), 과염소산 O (O_3 vs O_2),
  우라실류(고리 밖 C=O) 실패. 의도적 차이: S–H 의 H_HB(dreiding.lt 에 S 주개 라벨 없음),
  페놀 O(O_3 유지 — 논문 이면각 규칙 (i)가 O_3–X_R 을 전제), 티오우레아 S_2(dreiding.lt 에 없음).
- moltemplate 는 결합 차수를 라벨 쌍으로만 정하므로 C_R–O_2 는 1.0 으로 들어감(논문은 1.5).
  O_R 로 두면 차수 1.5 지만 r0 가 0.1 Å 길어짐(카복실레이트 C–O 1.35 vs 실제 1.25).
- dreiding.lt 파라미터: S D0 0.32(두 전사본은 0.344), Cl 질량 34.45, σ 산술평균 혼합,
  improper 는 중심당 1개(K=40), 수소결합 D_hb 4.0(Gasteiger 전하와 맞춘 값)으로 확인.
  결합 r0 는 Fe 를 빼고 전부 R_I+R_J−0.01 과 일치, 평형각도 In 외 일치.
- 조용히 틀리는 경우: ClO4⁻ 의 Cl 중심 평형각 180°(Cl 은 1가 전용 타입).
  다이아릴아이오도늄은 이면각 규칙이 없어 오류로 멈춤(메시지 불친절).

## 5단계 — 2차 검토 (2026-10-31)
### 1차 검토 정정
- 수소결합 D_hb 4.0 을 "Gasteiger 전하 기준 값"이라고 한 것은 근거 없음. 논문을 인용한 LUNAR 전사본은
  Gasteiger 7.0 / 전하 없음 9.0 kcal/mol, Goddard 강의 자료도 Gasteiger 7.0, LAMMPS 문서 예시는 9.5.
  moltemplate 생성 스크립트(Dreiding_nonbondHB.py)는 4.0 의 출처를 밝히지 않음.
- 술폰아마이드 N → N_3 은 확실한 수정이 아니라 절충: N 각도는 좋아지지만(+9.9° → −3.1°) S–N 길이가
  나빠짐(+0.005 → +0.057 Å, MMFF94 대비).
- 아마이드 C_R/N_R 은 C–N(1.320 → 1.349, MMFF 1.366)·장벽(5 → 25, 실험 16~21)은 좋아지나 C=O 가
  1.222 → 1.255 로 길어짐(MMFF 1.223). 또 양쪽에 방향족이 붙은 아마이드(벤즈아닐라이드, NHPI,
  폴리이미드·아라미드 단위)는 b1 플래그 충돌로 라벨 조합이 없음 → 그룹별 되돌림(fallback) 필요.
- dreid-typer 는 정답지가 아님: 술포네이트 O 타입이 SMILES 표기(어느 O 가 O⁻ 인지)에 따라 바뀜,
  P=O 를 O_3 으로 둠, TCM⁻ 중심 C 를 C_3 으로 둠.
### 근거 보강
- S ε 0.344: SEAMM, LUNAR, lammps_interface, pysimm 네 전사본 모두 일치 (moltemplate 만 0.32).
- 혼합 규칙: LUNAR(@combination geometric), Towhee 모두 기하평균. moltemplate 는 σ 산술평균.
- improper: lammps_interface 는 평면 중심당 3개 × K/3, moltemplate 는 1개 × K=40 → 대칭 중심에서 동등.
### 새로 찾은 것
- cleanup_moltemplate.sh 를 쓰면 수소결합 주개 플래그가 틀어짐: 정리 후 타입 번호가 받개 < 주개가 되면
  `pair_coeff 받개 주개 ... i` 로 남아 LAMMPS 가 받개를 주개로 봄(H 가 없으니 항이 0).
  메탄올+아세톤 MD 에서 수소결합 에너지 −77.4 → −3.8 kcal/mol. 정리 전(원래 출력)은 정상.
  moltemplate postprocess 는 pair_style 을 알 때만 i↔j 를 바꾸는데, 정리 경로에서 그게 안 되는 것으로 보임.
- 할라이드 음이온(Cl⁻, Br⁻, I⁻)은 공유결합 할로젠 LJ(σ 3.5~3.7 Å)를 받음. 이온용 값(σ 4.4~5.4)보다 작음.
- 생성 .lt 에서 보정 가능 여부: `pair_coeff` 덮어쓰기(S ε)는 정리 전·후 모두 유지,
  `mass` 덮어쓰기(Cl)는 정리 전만 유지(정리 스크립트가 지움).
- 84개 분자(현재·제안 규칙)를 moltemplate 로 돌려 각·이면각·improper 누락 0 확인.
  MMFF94 대비 DREIDING 고유 오차: S–O +0.13 Å, TFSI S–N +0.10, C(=O)–C(sp3) −0.06,
  카복실레이트 O–C–O −10°, 술포닐 각 ±6~10°.
- RDKit MMFFOptimizeMolecule 은 분자의 방향족 표시를 바꿈 → 타이핑 전에 쓰면 안 됨
  (CLI 는 조각 복사본에만 써서 영향 없음을 확인).
- 제안 규칙 시제품(스크래치)으로 dreid-typer 와의 차이 89 → 52 원자(188 분자), 남은 차이는 의도적.
