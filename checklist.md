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

## 2단계 (2026-10-26)
- [x] 매뉴얼 예시 라벨과 도구 라벨이 moltemplate 결과(모든 결합·각·이면각 계수)에서 같은지 테스트
- [취소] moltemplate.sh 실행·cleanup·LAMMPS 입력 보정 — 사용자는 .lt 까지만 원함

## 3단계 — 오프라인 HTML 도구 (2026-10-28)
- [x] dreiding.lt 규칙표 추출 (web/build_tables.py)
- [x] 타이퍼 JS 이식 + Python 대조 테스트 (Node)
- [x] mol2/pdb/xyz 파서, SMILES 그래프 대응, .lt 작성 (Python CLI 출력과 동일성 테스트)
- [x] SMILES만일 때 JS 3D 생성 (DREIDING 결합·각 + 반발 최소화)
- [x] 화면: 입력, 경고, 타입 요약표, 3D 뷰어(라벨·호버·강조·색 모드), .lt 출력(복사·다운로드)
- [x] RDKit.js wasm 압축 내장 단일 HTML 빌드 (web/build.py)
- [x] Playwright 로 file:// 오프라인 로딩·입력별 결과 확인
- [x] README·메모 갱신, 커밋·push
- [x] Playwright 브라우저 테스트를 저장소에 추가 (web/tests/test_browser.js)
