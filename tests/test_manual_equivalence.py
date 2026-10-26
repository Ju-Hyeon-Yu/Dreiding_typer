# 라벨 매뉴얼 예시 라벨과 도구 라벨이 moltemplate 결과 계수에서 같은지 확인하는 테스트
import os
import re
import shutil
import subprocess
import sys

import pytest
from rdkit import Chem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import dreiding_lt  # noqa: E402
from dreiding_typer import DreidingFF, assign_labels  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("moltemplate.sh") is None, reason="moltemplate 없음")

# DREIDING_Label_Manual.pdf 그림 1b, 2a 의 라벨 (무거운 원자 SMILES 순서, H 는 규칙대로)
MANUAL = {
    "cpd": ("O=C1C=CC=C1", ["O_2_b2", "C_2_b2", "C_2", "C_2_b1", "C_2_b1", "C_2"]),
    "imid": ("c1cnc[nH]1", ["C_R", "C_R", "N_R_d1_ha", "C_R", "N_R_b1_d2_hd"]),
}


def interactions(d, n_by_mol):
    """분자별로 (상호작용 종류, 원자 인덱스) → 계수 문자열."""
    data = open(d / "system.data").read()
    settings = [l.split() for l in open(d / "system.in.settings") if l.strip()]
    coeff = {}
    for t in settings:
        if t[0] in ("bond_coeff", "angle_coeff", "dihedral_coeff", "improper_coeff"):
            coeff[(t[0].split("_")[0], int(t[1]))] = " ".join(t[2:])
    masses = {int(m.group(1)): float(m.group(2))
              for m in re.finditer(r"^\s*(\d+)\s+([\d.]+)", data.split("Masses")[1].split("Atoms")[0], re.M)}
    out, start = [], 0
    bounds = []
    for n in n_by_mol:
        bounds.append((start + 1, start + n))
        start += n

    def owner(aid):
        for k, (lo, hi) in enumerate(bounds):
            if lo <= aid <= hi:
                return k, aid - lo
        raise KeyError(aid)

    res = [dict() for _ in n_by_mol]
    atoms = re.search(r"\nAtoms[^\n]*\n\n(.*?)\n\n", data, re.S).group(1).splitlines()
    for line in atoms:
        t = line.split()
        k, i = owner(int(t[0]))
        res[k][("mass", i)] = masses[int(t[2])]
    for sec, kind, n in (("Bonds", "bond", 2), ("Angles", "angle", 3),
                         ("Dihedrals", "dihedral", 4), ("Impropers", "improper", 4)):
        blk = re.search(rf"\n{sec}[^\n]*\n\n(.*?)(\n\n|$)", data, re.S)
        for line in blk.group(1).splitlines():
            t = line.split()
            ids = [owner(int(x)) for x in t[2:2 + n]]
            k = ids[0][0]
            idx = tuple(i for _, i in ids)
            if kind == "improper":
                idx = (idx[0], tuple(sorted(idx[1:])))
            elif idx[0] > idx[-1]:
                idx = idx[::-1]
            res[k][(kind, idx)] = coeff[(kind, int(t[1]))]
    return res


@pytest.mark.parametrize("key", list(MANUAL))
def test_manual_labels_give_same_parameters(key, tmp_path):
    smiles, heavy = MANUAL[key]
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ours, _ = assign_labels(mol, DreidingFF())
    manual = list(ours)
    for a, lab in zip([a for a in mol.GetAtoms() if a.GetSymbol() != "H"], heavy):
        manual[a.GetIdx()] = lab
    assert ours != manual, "매뉴얼과 라벨 이름이 같으면 이 테스트는 의미가 없음"
    xyz = dreiding_lt.embed(mol, seed=1)
    bonds = sorted(tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds())
    for tag, labs in (("ours", ours), ("man", manual)):
        rows = [{"name": f"{a.GetSymbol()}{a.GetIdx() + 1}", "frag": 1, "label": labs[a.GetIdx()],
                 "q": 0.0, "xyz": xyz[a.GetIdx()]} for a in mol.GetAtoms()]
        dreiding_lt.write_lt(tmp_path / f"{tag}.lt", tag, smiles, rows, bonds, "0")
    (tmp_path / "system.lt").write_text(
        'import "ours.lt"\nimport "man.lt"\na = new ours\nb = new man.move(30,0,0)\n'
        'write_once("Data Boundary") {\n -50 80 xlo xhi\n -50 50 ylo yhi\n -50 50 zlo zhi\n}\n')
    r = subprocess.run(["moltemplate.sh", "system.lt"], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-1500:]
    a, b = interactions(tmp_path, [mol.GetNumAtoms()] * 2)
    assert a.keys() == b.keys()
    diff = {k: (a[k], b[k]) for k in a if a[k] != b[k]}
    assert not diff, diff
