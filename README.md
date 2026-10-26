# Dreiding_typer

SMILES(필요하면 RESP 전하가 든 mol2와 함께)에서 moltemplate DREIDING 단일 분자 `.lt` 파일을
만드는 도구입니다. 원자 타입은 moltemplate 배포본 `force_fields/dreiding.lt`(Matthew Bone 판)의
라벨 체계로 자동 할당하고, 결과 파일은 `import "dreiding.lt"` + `MOL inherits DREIDING { ... }`
형태라 `system.lt`에서 바로 불러 쓸 수 있습니다.

- 타입은 SMILES의 결합 차수로 정하므로 구조(conformer)에 따라 바뀌지 않습니다.
- 수소결합 항(`H_HB`, `_hd`, `_ha`)을 자동으로 표시합니다.
- 결합 차수 플래그(`_b1`/`_b2`)는 `dreiding.lt`의 규칙을 직접 읽어, moltemplate가 붙일
  결합 차수가 실제와 같고 모든 결합에 이면각 규칙이 있는 조합으로 고릅니다.

## 설치

```bash
git clone https://github.com/Ju-Hyeon-Yu/Dreiding_typer.git
cd Dreiding_typer
pip install rdkit moltemplate
```

Python 3.9 이상, RDKit, moltemplate가 필요합니다. 라벨 목록과 규칙은 설치된 moltemplate의
`dreiding.lt`에서 실행할 때마다 읽습니다.

## 사용법

```bash
# RESP 전하가 든 mol2 사용 — 원자 순서·좌표·전하는 mol2 그대로, 결합 차수는 SMILES에서
python3 dreiding_lt.py --smiles "Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F" \
    --mol2 ip.mol2 -o single_molecule.lt --name MOL

# SMILES만 — RDKit 3D 구조, 전하 0 (경고 표시)
python3 dreiding_lt.py --smiles "CCO" -o ethanol.lt --name MOL

# 여러 분자는 셸 반복문으로 (list.csv: 이름,SMILES)
while IFS=, read name smi; do
  python3 dreiding_lt.py --smiles "$smi" --mol2 "$name.mol2" -o "$name.lt" --name "$name"
done < list.csv
```

| 옵션 | 설명 |
|---|---|
| `--smiles` | 필수. 이온쌍은 `.`으로 연결 (예: `A.B`) |
| `--mol2` | 좌표·전하를 가져올 mol2. 없으면 RDKit으로 3D 구조를 만들고 전하는 0 |
| `-o`, `--out` | 출력 `.lt` 경로 |
| `--name` | moltemplate 객체 이름. 기본값은 출력 파일 이름 |
| `--dreiding-lt` | 라벨·규칙을 읽을 `dreiding.lt` 경로. 기본값은 설치된 moltemplate의 것 |
| `--seed` | mol2가 없을 때 3D 구조 생성 난수 시드 |

화면에는 원자별 라벨과 전하가 출력되고, DREIDING으로 정확히 표현할 수 없는 부분은 `경고:`로
표준 오류에 나옵니다. 라벨을 정할 수 없으면 `.lt`를 쓰지 않고 멈춥니다.

### 출력 예

```
# dreiding_lt.py 생성 — SMILES: Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F
# 전하: mol2 (ip.mol2)

import "dreiding.lt"

MOL inherits DREIDING {

  write("Data Atoms") {
    $atom:C1     $mol:f1 @atom:C_3              0.079135     3.44930     1.07820     0.99060
    $atom:N2     $mol:f1 @atom:N_R_d2          -0.239526     2.56160    -0.01380     0.61200
    ...
  }

  write("Data Bond List") {
    $bond:b1    $atom:C1     $atom:N2
    ...
  }

}
```

이온쌍은 조각마다 다른 molecule ID(`$mol:f1`, `$mol:f2`)가 붙습니다.

### system.lt에서 불러오기

```
import "single_molecule.lt"
m = new MOL
```

이후 `moltemplate.sh system.lt`로 LAMMPS 입력을 만들면 됩니다. `dreiding.lt` 기본 설정은
`hybrid/overlay hbond/dreiding/lj ... lj/cut/coul/long` + `kspace_style ewald`이고,
LAMMPS MOLECULE 패키지(`hbond/dreiding/lj`, `dihedral charmm`, `improper umbrella`)가 필요합니다.

## SMILES가 없고 mol2·PDB만 있을 때

`--smiles`는 필수입니다. DREIDING 라벨(`C_R`/`C_2`, `_b1` 등)이 결합 차수로 정해지는데,
결합 차수와 형식 전하를 SMILES에서만 가져오기 때문입니다(mol2의 결합 차수 열은 쓰지 않음).

SMILES가 없으면 아래처럼 좌표에서 결합 차수를 추정해 SMILES를 먼저 만든 뒤 넣으면 됩니다.
RDKit `DetermineBonds`를 쓰며, 파일의 결합 정보는 무시하고 좌표만 사용합니다.

```python
# guess_smiles.py — 사용: python3 guess_smiles.py ip.mol2 0   (마지막 인자는 계 전체 총 전하)
import sys
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

path, charge = sys.argv[1], int(sys.argv[2])
if path.endswith(".pdb"):
    mol = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False, proximityBonding=False)
else:
    mol = Chem.MolFromMol2File(path, removeHs=False, sanitize=False, cleanupSubstructures=False)
rw = Chem.RWMol(mol)
for b in list(rw.GetBonds()):
    rw.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
for a in rw.GetAtoms():
    a.SetFormalCharge(0)
    a.SetNoImplicit(True)
mol = rw.GetMol()
rdDetermineBonds.DetermineBonds(mol, charge=charge)
print(Chem.MolToSmiles(Chem.RemoveHs(mol)))
```

```bash
smi=$(python3 guess_smiles.py ip.mol2 0)
echo "$smi"                                   # 반드시 눈으로 확인
python3 dreiding_lt.py --smiles "$smi" --mol2 ip.mol2 -o single_molecule.lt --name MOL
```

주의할 점은 다음과 같습니다.

- **H가 모두 들어 있는 구조**여야 합니다. 총 전하는 직접 넣어야 합니다(이온쌍이면 보통 0).
- **공명 구조가 여러 개인 이온은 다른 표기로 나올 수 있습니다.** 예를 들어 TFSI⁻는 음전하가
  N이 아닌 O에, S=N 이중결합으로 추정되어 N 라벨이 `N_R_d1` 대신 `N_2_d1`이 됩니다. 출력된
  SMILES가 의도한 구조와 다르면 직접 고쳐서 넣으십시오.
- **좌표와 전하는 mol2에서만 읽습니다.** PDB만 있으면 위 방법으로 SMILES는 얻을 수 있지만,
  `--mol2` 없이 실행하면 좌표는 RDKit이 새로 만들고 전하는 0이 됩니다. PDB 좌표와 RESP 전하를
  쓰려면 mol2로 바꿔서 넣으십시오.

## mol2 요구 사항

- 원소는 원자 타입 열(SYBYL `C.ar`, GAFF `ca` 등)에서 읽고, 안 되면 원자 이름 앞 글자로 읽습니다.
- 전하는 9번째 열을 씁니다. 없으면 0으로 두고 경고합니다.
- **원자 순서는 맞출 필요가 없습니다.** SMILES와 mol2의 원자를 원소+연결 구조(그래프 동형)로
  자동 대응시키고, 출력 `.lt`는 mol2의 원자 순서를 그대로 따릅니다. 그래서 mol2의 결합 차수
  열은 쓰지 않고 BOND 섹션의 연결 정보만 씁니다. 원자·결합 수나 연결 구조가 다르면 멈춥니다.
- 전하 합이 SMILES의 형식 전하와 다르면 경고합니다.

## 타입 규칙 요약

| 대상 | 라벨 |
|---|---|
| 방향족 C/N/O | `C_R`, `N_R_d1`(이웃 2), `N_R_d2`(이웃 3), `O_R` |
| 카복실레이트·니트로 공명 묶음 | 중심 `C_R`/`N_R`, O는 `O_R` (결합 차수 1.5) |
| sp2/방향족 원자에 붙은 3배위 N (아마이드, 아닐린) | `N_R_d2` |
| 다중결합 없는 2배위 N (이미드 음이온) | `N_R_d1` |
| 에스터·페놀·산 OH의 O | `O_3` (원 논문 이면각 규칙 (i) 적용) |
| S | `S_3` (유일한 S 타입) |
| I, Br, Cl, F | `I`, `Br`, `Cl`, `F` |
| N/O/F에 붙은 H | `H_HB`, 그 무거운 원자는 `_hd` |
| 수소결합 받개 | O, F, 비공유전자쌍이 남는 N에 `_ha` |
| 단일결합으로 이어진 sp2 원자 쌍 | `_b1`/`_b2`를 탐색으로 부여 |

### 라벨 매뉴얼 예시와의 관계

모든 라벨은 moltemplate `dreiding.lt`에 정의된 것만 씁니다. 매뉴얼
([DREIDING_Label_Manual.pdf](https://github.com/m-bone/moltemplate-DREIDING)) 그림 예시와
이름이 다른 경우가 있지만, moltemplate가 붙이는 결합·각·이면각·improper 계수가 모두 같음을
테스트로 확인했습니다(사이클로펜타디에논, 이미다졸 N–H).

예외는 아마이드 하나입니다. 매뉴얼 예시(`C_2` + `N_2_b1_d2_hd`)는 C–N을 이중결합
(평형 1.275 Å, 회전 장벽 45 kcal/mol)으로 넣고, 이 도구(`C_2_b1` + `N_R_b1_d2_hd`)는
단일결합(1.31 Å, 장벽 5 kcal/mol)으로 넣습니다. 실제 아마이드(약 1.33 Å, 15–20 kcal/mol)에
더 가까운 쪽을 택했습니다.

## 알려진 한계 (실행 시 경고로 표시)

- **3·4배위 S(술포네이트, 술폰, TFSI)** — DREIDING에 해당 타입이 없어 `S_3`(H₂S 기준,
  평형각 92.1°)를 씁니다. S=O는 단일결합 힘상수로 들어가고, S를 지나는 이면각 장벽은
  이웃 2개를 가정해 나누므로 실제보다 3배 큽니다.
- **4배위 N(4급 암모늄, 양성자화 아민)** — `N_3`가 이웃 3개를 가정해 이 N을 지나는
  이면각 장벽이 1.5배로 들어갑니다.
- **Cl** — `dreiding.lt`의 Cl 질량이 34.45로 들어가 있습니다(실제 35.45).

## 파일 구성

| 파일 | 내용 |
|---|---|
| `dreiding_lt.py` | CLI — SMILES/mol2 입력, 원자 대응, `.lt` 쓰기 |
| `dreiding_typer.py` | 타이퍼 — 라벨 규칙, `dreiding.lt` 규칙 파서, `_b1`/`_b2` 탐색 |
| `tests/` | 단위·통합 테스트 |
| `plan.md`, `checklist.md`, `context-notes.md` | 설계 계획, 작업 목록, 결정과 근거 기록 |

## 테스트

```bash
pip install pytest lammps        # lammps 는 LAMMPS 단일점 테스트용 (없으면 그 테스트만 건너뜀)
python3 -m pytest -q tests       # moltemplate 실행 포함 약 1분
```

통합 테스트는 이온쌍, TFSI/BMIM, 아마이드, 카복실레이트, 니트로, 비페닐, 공액 사슬 등
14개 분자의 `.lt`를 moltemplate로 돌려 다음을 확인합니다.

- 모든 결합의 힘상수가 350 × 기대 결합 차수와 같음
- 모든 결합각·이면각이 타입을 받음 (누락 없음)
- 이면각 계수가 Mayo 1990 규칙 (a)~(j)의 독립 구현과 같음 (위 한계 항목 제외)
- 이온쌍 조각마다 molecule ID가 다르고, mol2 원자 순서·전하가 유지됨
- 매뉴얼 예시 라벨과 도구 라벨의 moltemplate 계수가 같음
- LAMMPS가 결과를 읽고 단일점 계산을 할 수 있음

## 참고

- Mayo, Olafson, Goddard, *J. Phys. Chem.* **94**, 8897 (1990) — DREIDING
- [jewettaij/moltemplate](https://github.com/jewettaij/moltemplate) — `force_fields/dreiding.lt`
- [m-bone/moltemplate-DREIDING](https://github.com/m-bone/moltemplate-DREIDING) — 라벨 매뉴얼
