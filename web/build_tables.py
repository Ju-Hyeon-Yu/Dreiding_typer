# dreiding.lt 에서 브라우저용 DREIDING 규칙표(라벨·결합 차수·r0·이면각·평형각)를 JSON 으로 추출
"""사용: python3 web/build_tables.py [dreiding.lt 경로] > web/build/tables.json

라벨 쌍별 결합 차수와 이면각 규칙 존재 여부는 dreiding_typer.DreidingFF 로 계산하므로
Python 타이퍼와 판정이 같다. r0 와 평형각은 3D 좌표 생성(최소화)에만 쓴다.
"""
import fnmatch
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dreiding_typer import DreidingFF  # noqa: E402


def read_sections(path):
    """write_once 블록별 줄 목록."""
    secs, cur = {}, None
    with open(path) as f:
        for line in f:
            s = line.strip()
            m = re.match(r'write_once\("([^"]+)"\)', s)
            if m:
                cur = m.group(1)
                secs.setdefault(cur, [])
                continue
            if s.startswith("}"):
                cur = None
                continue
            if cur and s and not s.startswith("#"):
                secs[cur].append(s.split())
    return secs


def last_match(rules, names):
    """moltemplate By Type 규칙: 나중 규칙이 이긴다."""
    hit = None
    for rule_name, pats in rules:
        if all(fnmatch.fnmatchcase(n, p) for n, p in zip(names, pats)) or \
                all(fnmatch.fnmatchcase(n, p) for n, p in zip(names[::-1], pats)):
            hit = rule_name
    return hit


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else None
    ff = DreidingFF(path)
    secs = read_sections(ff.path)
    settings = [t for k, v in secs.items() if k == "In Settings" for t in v]
    bond_r0 = {t[1].split(":", 1)[1]: float(t[3]) for t in settings if t[0] == "bond_coeff"}
    angle_t0 = {t[1].split(":", 1)[1]: float(t[3]) for t in settings if t[0] == "angle_coeff"}
    bond_rules = [(t[0].split(":", 1)[1], [x.split(":", 1)[1] for x in t[1:3]])
                  for t in secs["Data Bonds By Type"]]
    angle_rules = [(t[0].split(":", 1)[1], [x.split(":", 1)[1] for x in t[1:4]])
                   for t in secs["Data Angles By Type"]]

    labels = sorted(ff.labels)
    pairs = {}
    for i, a in enumerate(labels):
        for b in labels[i:]:
            order = ff.bond_order(a, b)
            if order is None:
                continue
            r0 = bond_r0[last_match(bond_rules, [a, b])]
            pairs[f"{a}|{b}"] = [round(order, 3), r0, 1 if ff.has_dihedral(a, b) else 0]
    theta0 = {}
    for c in labels:
        rule = last_match(angle_rules, ["X", c, "X"])
        if rule:
            theta0[c] = angle_t0[rule]
    out = {
        "source": os.path.basename(ff.path),
        "labels": {k: ff.labels[k] for k in labels},
        "pairs": pairs,          # "A|B" (A<=B) → [결합 차수, r0, 이면각 규칙 유무]
        "theta0": theta0,        # 중심 라벨 → 평형각(도)
    }
    json.dump(out, sys.stdout, separators=(",", ":"), ensure_ascii=False)


if __name__ == "__main__":
    main()
