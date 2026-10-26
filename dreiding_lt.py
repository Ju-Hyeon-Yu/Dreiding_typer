#!/usr/bin/env python3
# SMILES(+RESP mol2)에서 moltemplate DREIDING .lt 파일을 만드는 CLI
"""사용 예
  python3 dreiding_lt.py --smiles "Cn1cc[n+](C)c1.[O-]S(=O)(=O)C(F)(F)F" --mol2 ip.mol2 -o ip.lt
  python3 dreiding_lt.py --smiles "CCO" -o ethanol.lt          # 좌표는 RDKit, 전하 0

mol2 를 주면 원자 순서·좌표·전하는 mol2 그대로 쓰고, 결합 차수와 화학 정보는
SMILES 에서 가져온다(두 분자를 원소+연결 그래프로 대응).
"""
import argparse
import os
import re
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

from dreiding_typer import DreidingFF, TypingError, assign_labels, default_lt_path

ELEMENTS = {Chem.GetPeriodicTable().GetElementSymbol(z) for z in range(1, 119)}
GAFF_PREFIX = [("cl", "Cl"), ("br", "Br"), ("i", "I"), ("f", "F"), ("c", "C"), ("n", "N"),
               ("o", "O"), ("s", "S"), ("p", "P"), ("h", "H")]


def _element(name, atype):
    """mol2 원자 타입(SYBYL 또는 GAFF) → 원소, 안 되면 원자 이름의 앞 글자."""
    base = atype.split(".")[0]
    if base in ELEMENTS:
        return base
    if atype.islower() or (atype[:1].islower()):
        for pre, el in GAFF_PREFIX:
            if atype.startswith(pre):
                return el
    letters = re.match(r"[A-Za-z]+", name).group(0)
    for cand in (letters[:2].capitalize(), letters[:1].upper()):
        if cand in ELEMENTS:
            return cand
    raise ValueError(f"mol2 원자 '{name} {atype}' 의 원소를 알 수 없습니다")


def read_mol2(path):
    atoms, bonds, sec = [], [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                sec = s[9:]
                continue
            if not s or s.startswith("#"):
                continue
            t = s.split()
            if sec == "ATOM":
                charge = float(t[8]) if len(t) > 8 else 0.0
                atoms.append({"id": int(t[0]), "name": t[1], "xyz": tuple(map(float, t[2:5])),
                              "el": _element(t[1], t[5]), "q": charge, "has_q": len(t) > 8})
            elif sec == "BOND":
                bonds.append((int(t[1]), int(t[2])))
    if not atoms:
        raise ValueError(f"{path}: ATOM 섹션이 비어 있습니다")
    return atoms, bonds


def _graph_mol(elements, bonds):
    """원소와 단일결합만 가진 비교용 분자 (결합 차수·전하 무시)."""
    m = Chem.RWMol()
    for el in elements:
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        m.AddAtom(a)
    for i, j in bonds:
        m.AddBond(i, j, Chem.BondType.SINGLE)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def match_mol2(template, atoms, bonds):
    """SMILES 템플릿 원자 인덱스 → mol2 원자 순번(0부터) 대응."""
    idx = {a["id"]: k for k, a in enumerate(atoms)}
    tgt = _graph_mol([a["el"] for a in atoms], [(idx[i], idx[j]) for i, j in bonds])
    qry = _graph_mol([a.GetSymbol() for a in template.GetAtoms()],
                     [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in template.GetBonds()])
    if tgt.GetNumAtoms() != qry.GetNumAtoms() or tgt.GetNumBonds() != qry.GetNumBonds():
        raise ValueError(f"원자/결합 수가 다릅니다: SMILES {qry.GetNumAtoms()}/{qry.GetNumBonds()}, "
                         f"mol2 {tgt.GetNumAtoms()}/{tgt.GetNumBonds()}")
    match = tgt.GetSubstructMatch(qry)
    if not match:
        raise ValueError("SMILES 와 mol2 의 연결 구조가 일치하지 않습니다 (원소 판정이나 결합 확인 필요)")
    return list(match)


def embed(template, seed):
    m = Chem.Mol(template)
    if AllChem.EmbedMolecule(m, randomSeed=seed) != 0:
        raise ValueError("RDKit 3D 구조 생성에 실패했습니다")
    if AllChem.MMFFHasAllMoleculeParams(m):
        AllChem.MMFFOptimizeMolecule(m)
    else:
        AllChem.UFFOptimizeMolecule(m)
    return [tuple(m.GetConformer().GetAtomPosition(i)) for i in range(m.GetNumAtoms())]


def write_lt(path, name, smiles, rows, bonds, charge_note):
    with open(path, "w") as f:
        f.write(f"# dreiding_lt.py 생성 — SMILES: {smiles}\n# 전하: {charge_note}\n\n")
        f.write('import "dreiding.lt"\n\n')
        f.write(f"{name} inherits DREIDING {{\n\n")
        f.write('  write("Data Atoms") {\n')
        for r in rows:
            x, y, z = r["xyz"]
            f.write(f"    $atom:{r['name']:<6} $mol:f{r['frag']} @atom:{r['label']:<14} "
                    f"{r['q']:>10.6f} {x:>11.5f} {y:>11.5f} {z:>11.5f}\n")
        f.write("  }\n\n")
        f.write('  write("Data Bond List") {\n')
        for k, (i, j) in enumerate(bonds, 1):
            f.write(f"    $bond:b{k:<4} $atom:{rows[i]['name']:<6} $atom:{rows[j]['name']}\n")
        f.write("  }\n\n}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description="SMILES(+RESP mol2) → moltemplate DREIDING .lt")
    ap.add_argument("--smiles", required=True, help="이온쌍은 '.' 으로 연결 (예: A.B)")
    ap.add_argument("--mol2", help="좌표·전하를 가져올 mol2 (원자 순서 유지)")
    ap.add_argument("-o", "--out", required=True, help="출력 .lt 경로")
    ap.add_argument("--name", help="moltemplate 분자 이름 (기본: 출력 파일 이름)")
    ap.add_argument("--dreiding-lt", default=None, help=f"라벨 목록을 읽을 dreiding.lt (기본: {default_lt_path()})")
    ap.add_argument("--seed", type=int, default=1, help="mol2 없을 때 3D 구조 난수 시드")
    args = ap.parse_args(argv)

    name = args.name or re.sub(r"\W", "_", os.path.splitext(os.path.basename(args.out))[0])
    mol0 = Chem.MolFromSmiles(args.smiles)
    if mol0 is None:
        sys.exit("SMILES 를 읽지 못했습니다")
    template = Chem.AddHs(mol0)
    n = template.GetNumAtoms()
    warnings = []

    if args.mol2:
        atoms, mbonds = read_mol2(args.mol2)
        try:
            t2m = match_mol2(template, atoms, mbonds)
        except ValueError as e:
            sys.exit(f"mol2 대응 실패: {e}")
        m2t = {m: t for t, m in enumerate(t2m)}
        # 원자 번호를 mol2 순서로 바꿔서 이후 출력·경고 번호가 mol2 와 같게 함
        template = Chem.RenumberAtoms(template, [m2t[k] for k in range(n)])
        xyz = [a["xyz"] for a in atoms]
        q = [a["q"] for a in atoms]
        if not all(a["has_q"] for a in atoms):
            warnings.append("mol2 일부 원자에 전하 열이 없어 0 으로 둡니다")
        charge_note = f"mol2 ({os.path.basename(args.mol2)})"
        qsum, formal = sum(q), Chem.GetFormalCharge(template)
        if abs(qsum - formal) > 1e-3:
            warnings.append(f"mol2 전하 합 {qsum:.4f} 이 SMILES 형식 전하 {formal} 와 다릅니다")
    else:
        xyz = embed(template, args.seed)
        q = [0.0] * n
        charge_note = "없음 — 모두 0 (RESP 는 --mol2 로 넣으세요)"
        warnings.append("mol2 가 없어 전하를 모두 0 으로 썼습니다")

    try:
        labels, w = assign_labels(template, DreidingFF(args.dreiding_lt))
    except TypingError as e:
        sys.exit(f"타이핑 실패: {e}")
    warnings += w

    frags = Chem.GetMolFrags(template)
    frag_of = {i: k for k, fr in enumerate(frags, 1) for i in fr}
    rows = [{"name": f"{a.GetSymbol()}{a.GetIdx() + 1}", "frag": frag_of[a.GetIdx()],
             "label": labels[a.GetIdx()], "q": q[a.GetIdx()], "xyz": xyz[a.GetIdx()]}
            for a in template.GetAtoms()]
    bonds = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in template.GetBonds())
    write_lt(args.out, name, args.smiles, rows, bonds, charge_note)

    width = max(len(r["label"]) for r in rows)
    print(f"{args.out}: 원자 {n}개, 결합 {len(bonds)}개, 조각 {len(frags)}개")
    for k, r in enumerate(rows, 1):
        print(f"  {k:4d} {r['name']:<6} {r['label']:<{width}} {r['q']:>10.6f}")
    for w in dict.fromkeys(warnings):
        print(f"경고: {w}", file=sys.stderr)


if __name__ == "__main__":
    main()
