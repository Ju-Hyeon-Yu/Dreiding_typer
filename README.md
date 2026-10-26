# dreiding_lt — SMILES(+RESP mol2) → moltemplate DREIDING `.lt`

moltemplate 배포본 `force_fields/dreiding.lt`(Matthew Bone 판)의 라벨 체계로
원자 타입을 자동 할당하고 `.lt` 파일을 만듭니다. 타입은 SMILES의 결합 차수로
정하므로 conformer에 따라 바뀌지 않습니다.

## 설치

```bash
pip install rdkit moltemplate        # 테스트까지 돌리려면 pytest 추가
```

## 사용법

```bash
# RESP 전하가 든 mol2 사용 (원자 순서·좌표·전하는 mol2 그대로)
python3 dreiding_lt.py --smiles "Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F" --mol2 ip.mol2 -o ip.lt

# SMILES만 (RDKit 3D 구조, 전하 0)
python3 dreiding_lt.py --smiles "CCO" -o ethanol.lt

# 여러 분자는 셸 반복문으로
while IFS=, read name smi; do
  python3 dreiding_lt.py --smiles "$smi" --mol2 "$name.mol2" -o "$name.lt"
done < list.csv
```

이온쌍은 SMILES를 `.`으로 이어 주면 조각마다 다른 molecule ID(`$mol:f1`, `$mol:f2`)가
붙습니다. 화면에 원자별 라벨이 출력되고, DREIDING으로 정확히 표현할 수 없는 부분은
`경고:`로 표준 오류에 나옵니다.

### moltemplate로 LAMMPS 입력 만들기

```bash
cat > system.lt <<'EOF'
import "ip.lt"
m = new ip
write_once("Data Boundary") {
  -30 30 xlo xhi
  -30 30 ylo yhi
  -30 30 zlo zhi
}
EOF
moltemplate.sh system.lt        # system.data, system.in.init, system.in.settings 생성
```

`dreiding.lt` 기본 설정은 `hybrid/overlay hbond/dreiding/lj ... lj/cut/coul/long` +
`kspace_style ewald`입니다. 진공에서 직접합으로 계산하려면 `system.in.init`과
`system.in.settings`의 `lj/cut/coul/long`을 `lj/cut/coul/cut`으로 바꾸고
`kspace_style` 줄을 지우면 됩니다. LAMMPS는 MOLECULE 패키지(`hbond/dreiding/lj`,
`dihedral charmm`, `improper umbrella`)가 필요합니다.

## mol2 요구 사항

- 원소는 원자 타입 열(SYBYL `C.ar`, GAFF `ca` 등)에서 읽고, 안 되면 원자 이름 앞 글자로 읽습니다.
- 전하는 9번째 열을 씁니다. 없으면 0으로 두고 경고합니다.
- SMILES와 mol2는 원소+연결 구조로 대응시키므로 결합 차수 열은 무시합니다.
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
| N/O/F에 붙은 H | `H_HB`, 그 무거운 원자는 `_hd` |
| 수소결합 받개 | O, F, 비공유전자쌍이 남는 N에 `_ha` |
| 단일결합으로 이어진 sp2 원자 쌍 | `_b1`/`_b2`를 탐색으로 부여 |

`_b1`/`_b2`는 **모든 결합의 라벨 차수가 실제 차수와 같고, 모든 결합에 이면각 규칙이
있는** 조합 중 플래그가 가장 적은 것으로 고릅니다. 판정은 `dreiding.lt`의 규칙을 직접
읽어 하므로 moltemplate 결과와 같습니다. 만족하는 조합이 없으면 멈춥니다.

## 알려진 한계 (경고로 표시)

- **3·4배위 S(술포네이트, 술폰, TFSI)** — DREIDING에 해당 타입이 없어 `S_3`(H₂S 기준,
  평형각 92.1°)를 씁니다. S=O는 단일결합 힘상수로 들어가고, S를 지나는 이면각 장벽은
  이웃 2개를 가정해 나누므로 실제보다 3배 큽니다.
- **4배위 N(4급 암모늄, 양성자화 아민)** — `N_3`가 이웃 3개를 가정해 이 N을 지나는
  이면각 장벽이 1.5배로 들어갑니다.
- **Cl** — `dreiding.lt`의 Cl 질량이 34.45로 들어가 있습니다(실제 35.45).

## 테스트

```bash
python3 -m pytest -q tests      # 통합 테스트는 moltemplate 실행으로 약 20초
```

통합 테스트는 이온쌍, TFSI/BMIM, 아마이드, 카복실레이트, 니트로, 비페닐, 공액 사슬 등
14개 분자를 moltemplate로 돌려 다음을 확인합니다.

- 모든 결합의 힘상수가 350 × 기대 결합 차수와 같음
- 모든 결합각·이면각이 타입을 받음 (누락 없음)
- 이면각 계수가 Mayo 1990 규칙 (a)~(j)의 독립 구현과 같음 (위 한계 항목 제외)
- 이온쌍 조각마다 molecule ID가 다르고, mol2 원자 순서·전하가 유지됨
