# 체크리스트

- [x] dreiding.lt에서 허용 라벨 목록 읽기
- [x] 기본 타입 규칙 (원소·혼성·공명·d 플래그)
- [x] 수소결합 플래그 (H_HB, _hd, _ha)
- [x] 결합 차수 플래그 탐색 (_b1, _b2) + 불일치 시 오류
- [x] 결합 차수·이면각 규칙을 dreiding.lt에서 직접 읽어 판정 (DreidingFF)
- [x] 이면각 규칙이 없는 라벨 쌍 금지 (moltemplate가 조용히 빼는 문제)
- [x] 단위 테스트 — 매뉴얼 예제 분자
- [x] 단위 테스트 — 대상 작용기 분자 (이온, 아이오딘)
- [x] mol2 파서 + SMILES 그래프 대응
- [x] .lt 작성기 + CLI
- [x] SMILES만 줄 때 이온 조각 겹침 방지
- [x] 통합 테스트 — moltemplate 실행 후 결합 힘상수·이면각 누락·논문 규칙 대조
- [x] 통합 테스트 — LAMMPS 진공 단일점
- [x] 사용법 README

## 2단계 — moltemplate.sh 연결 (2026-10-26)
- [ ] `--lammps DIR`: system.lt 작성 → moltemplate.sh 실행 → cleanup_moltemplate.sh
- [ ] 수소결합 쌍이 없으면 pair_style 에서 hbond/dreiding/lj 제거 (LAMMPS 오류 방지)
- [ ] `--vacuum`: lj/cut/coul/cut, kspace 삭제, boundary s s s, 큰 cutoff
- [ ] 매뉴얼 예시 라벨과 도구 라벨이 moltemplate 결과(모든 결합·각·이면각 계수)에서 같은지 테스트
- [ ] LAMMPS 로 파이프라인 출력 단일점 확인 테스트 (주개 있는 계/없는 계, 진공/주기)
- [ ] README 갱신
