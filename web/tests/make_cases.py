# JS 대조용 입력(원자 순서를 섞은 mol2·pdb·xyz)과 Python CLI 의 .lt 정답을 만드는 스크립트
"""사용: python3 web/tests/make_cases.py OUTDIR

OUTDIR/cases.json 에 [{smiles, mol2, pdb, pdb_noconect, xyz, lt}] 경로 목록을 쓴다.
"""
import json
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from rdkit import Chem  # noqa: E402
import dreiding_lt  # noqa: E402

SMILES = [
    "Cn1cc[n+](Cc2ccc(I)cc2)c1.[O-]S(=O)(=O)C(F)(F)F",
    "FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F.CCCCn1cc[n+](C)c1",
    "CC(=O)[O-].C[NH3+]",
    "CC(=O)Nc1ccccc1",
    "O=C1C=CC=C1",
    "c1ccc(cc1)-c1ccccc1",
    "O=[N+]([O-])c1ccccc1",
    "C[N+](C)(C)C.[I-]",
    "CCCC[P+](CCCC)(CCCC)CCCC.[Br-]",
    "OC(=O)c1ccccc1C(=O)O",
    "Cc1ccc(cc1)S(=O)(=O)[O-].CCCC[N+]1(C)CCCC1",
    "N#C[N-]C#N.CCn1cc[n+](C)c1",
]


def write_files(smi, base, seed):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    xyz = dreiding_lt.embed(m, seed=seed)
    n = m.GetNumAtoms()
    order = list(range(n))
    random.Random(seed).shuffle(order)
    pos = {old: new for new, old in enumerate(order)}
    q = {}
    for new, old in enumerate(order):
        q[new] = round(0.001 * (old + 1) * (-1) ** old, 6)
    # mol2 (SYBYL 타입, 결합 차수는 모두 1 — 대응에는 연결만 씀)
    with open(base + ".mol2", "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\nt\n{n} {m.GetNumBonds()} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for new, old in enumerate(order):
            a = m.GetAtomWithIdx(old)
            x, y, z = xyz[old]
            sy = a.GetSymbol() + (".ar" if a.GetIsAromatic() else "")
            f.write(f"{new + 1} {a.GetSymbol()}{new + 1} {x:.4f} {y:.4f} {z:.4f} {sy} 1 MOL {q[new]:.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for k, b in enumerate(m.GetBonds(), 1):
            f.write(f"{k} {pos[b.GetBeginAtomIdx()] + 1} {pos[b.GetEndAtomIdx()] + 1} 1\n")
    # pdb / xyz (같은 원자 순서)
    pm = Chem.RenumberAtoms(m, order)
    conf = Chem.Conformer(n)
    for new, old in enumerate(order):
        conf.SetAtomPosition(new, xyz[old])
    pm.RemoveAllConformers()
    pm.AddConformer(conf)
    Chem.MolToPDBFile(pm, base + ".pdb")
    with open(base + ".pdb") as f:
        lines = [l for l in f if not l.startswith("CONECT")]
    with open(base + "_noconect.pdb", "w") as f:
        f.writelines(lines)
    with open(base + ".xyz", "w") as f:
        f.write(f"{n}\n{smi}\n")
        for new, old in enumerate(order):
            x, y, z = xyz[old]
            f.write(f"{m.GetAtomWithIdx(old).GetSymbol()} {x:.4f} {y:.4f} {z:.4f}\n")
    dreiding_lt.main(["--smiles", smi, "--mol2", base + ".mol2", "-o", base + ".lt", "--name", "MOL"])


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    cases = []
    saved = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = open(os.devnull, "w")
    try:
        for k, smi in enumerate(SMILES):
            base = os.path.join(out, f"case{k}")
            write_files(smi, base, seed=k + 3)
            cases.append({"smiles": smi, "mol2": base + ".mol2", "pdb": base + ".pdb",
                          "pdb_noconect": base + "_noconect.pdb", "xyz": base + ".xyz", "lt": base + ".lt"})
    finally:
        sys.stdout.close()
        sys.stdout, sys.stderr = saved
    with open(os.path.join(out, "cases.json"), "w") as f:
        json.dump(cases, f, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
