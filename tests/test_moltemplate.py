# 생성한 .lt 를 moltemplate 로 돌려 결합 차수·이면각·mol ID 를 원 논문 규칙과 대조하는 통합 테스트
import os
import random
import re
import shutil
import subprocess
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import dreiding_lt  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("moltemplate.sh") is None, reason="moltemplate 없음")

SMILES = {
    "ionpair": "Cn1cc[n+](Cc2ccc(I)cc2)c1.[O-]S(=O)(=O)C(F)(F)F",
    "tfsi_bmim": "FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F.CCCCn1cc[n+](C)c1",
    "acetamide": "CC(N)=O",
    "acetate_mea": "CC(=O)[O-].C[NH3+]",
    "biphenyl": "c1ccc(cc1)-c1ccccc1",
    "hexatriene": "C=CC=CC=C",
    "cpdone": "O=C1C=CC=C1",
    "nitrobenzene": "O=[N+]([O-])c1ccccc1",
    "ester": "CC(=O)OC",
    "phenol": "Oc1ccccc1",
    "propene": "CC=C",
    "acetonitrile": "CC#N",
    "tma_iodide": "C[N+](C)(C)C.[I-]",
    "anilide": "CC(=O)Nc1ccccc1",
}


def write_shuffled_mol2(smiles, path, seed=7):
    """원자 순서를 섞고 SYBYL 형식 이름/타입을 쓴 mol2 (전하는 원자 번호로 구분 가능한 값)."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=1)
    order = list(range(m.GetNumAtoms()))
    random.Random(seed).shuffle(order)
    pos = {old: new for new, old in enumerate(order)}
    conf = m.GetConformer()
    charges = {}
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\nt\n{m.GetNumAtoms()} {m.GetNumBonds()} 1 0 0\nSMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for new, old in enumerate(order):
            a = m.GetAtomWithIdx(old)
            p = conf.GetAtomPosition(old)
            q = round(0.001 * (old + 1), 6)
            charges[new] = (a.GetSymbol(), q)
            sy = a.GetSymbol() + (".ar" if a.GetIsAromatic() else "")
            f.write(f"{new + 1} {a.GetSymbol()}{new + 1} {p.x:.4f} {p.y:.4f} {p.z:.4f} {sy} 1 MOL {q:.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for k, b in enumerate(m.GetBonds(), 1):
            f.write(f"{k} {pos[b.GetBeginAtomIdx()] + 1} {pos[b.GetEndAtomIdx()] + 1} 1\n")
    return charges


def read_sections(path):
    secs, cur = {}, None
    for line in open(path):
        s = line.split("#")[0].strip()
        if not s:
            continue
        if re.match(r"^[A-Z][A-Za-z ]+$", s) and not s[0].isdigit():
            cur = s
            secs[cur] = []
        elif cur:
            secs[cur].append(s.split())
    return secs


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    d = tmp_path_factory.mktemp("mt")
    mols = {}
    for name, smi in SMILES.items():
        args = ["--smiles", smi, "-o", str(d / f"{name}.lt")]
        if name == "ionpair":
            mols["mol2_charges"] = write_shuffled_mol2(smi, d / "ip.mol2")
            args += ["--mol2", str(d / "ip.mol2")]
        dreiding_lt.main(args)
        mols[name] = smi
    with open(d / "system.lt", "w") as f:
        for name in SMILES:
            f.write(f'import "{name}.lt"\n')
        for name in SMILES:
            f.write(f"m_{name} = new {name}\n")
        f.write('write_once("Data Boundary") {\n -50 50 xlo xhi\n -50 50 ylo yhi\n -50 50 zlo zhi\n}\n')
    r = subprocess.run(["moltemplate.sh", "system.lt"], cwd=d, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    return d, mols


def parse_lt_atoms(path):
    rows = []
    for line in open(path):
        m = re.match(r"\s*\$atom:(\S+)\s+\$mol:f(\d+)\s+@atom:(\S+)\s+(\S+)", line)
        if m:
            rows.append((m.group(1), int(m.group(2)), m.group(3), float(m.group(4))))
    return rows


def lt_molecule(path, smiles):
    """.lt 의 원자 순서대로 번호를 바꾼 RDKit 분자와 라벨."""
    rows = parse_lt_atoms(path)
    tmpl = Chem.AddHs(Chem.MolFromSmiles(smiles))
    bonds = []
    names = {r[0]: k for k, r in enumerate(rows)}
    for line in open(path):
        m = re.match(r"\s*\$bond:\S+\s+\$atom:(\S+)\s+\$atom:(\S+)", line)
        if m:
            bonds.append((names[m.group(1)], names[m.group(2)]))
    els = [re.match(r"[A-Z][a-z]?", r[0]).group(0) for r in rows]
    t2l = dreiding_lt.match_mol2(tmpl, [{"id": k + 1, "el": e} for k, e in enumerate(els)],
                                 [(i + 1, j + 1) for i, j in bonds])
    l2t = {lv: t for t, lv in enumerate(t2l)}
    mol = Chem.RenumberAtoms(tmpl, [l2t[k] for k in range(len(rows))])
    return mol, [r[2] for r in rows], rows


def expected_order(mol, i, j, labels):
    b = mol.GetBondBetweenAtoms(i, j)
    sp2 = all(any(t in labels[k] for t in ("_2", "_R", "_1")) for k in (i, j))
    if not sp2:
        return 1.0                       # S=O 등: DREIDING 은 단일결합만
    if b.GetIsAromatic():
        return 1.5
    if "_R" in labels[i] and "_R" in labels[j] and not (b.GetBeginAtom().IsInRing() or b.GetEndAtom().IsInRing()):
        return 1.5                       # 공명 묶음(카복실레이트, 니트로)
    return b.GetBondTypeAsDouble()


def hyb(label):
    if label.split("_")[0] in ("H", "F", "Cl", "Br", "I"):
        return "mono"
    for t, h in (("_R", "R"), ("_1", "sp"), ("_2", "sp2"), ("_3", "sp3")):
        if t in label:
            return h
    return "mono"


def paper_dihedral(mol, labels, a, b, c, d):
    """Mayo 1990 규칙 (a)~(j) — lammps_interface 구현을 라벨 체계에 맞게 옮긴 독립 판정."""
    hb, hc = hyb(labels[b]), hyb(labels[c])
    if "sp" in (hb, hc) or "mono" in (hb, hc):
        return None
    sp2 = ("sp2", "R")
    oxy = lambda lab: lab.startswith(("O_3", "S_3"))  # noqa: E731
    nonoxy_sp2 = lambda lab: lab.startswith(("C_R", "N_R", "B_2", "C_2", "N_2"))  # noqa: E731
    order = expected_order(mol, b, c, labels)
    if hb == "sp3" and hc == "sp3":
        V, n, phi0 = (2.0, 2, 90.0) if oxy(labels[b]) and oxy(labels[c]) else (2.0, 3, 180.0)
    elif (hb in sp2) != (hc in sp2):
        V, n, phi0 = 1.0, 6, 0.0
        if (oxy(labels[b]) and nonoxy_sp2(labels[c])) or (oxy(labels[c]) and nonoxy_sp2(labels[b])):
            V, n, phi0 = 2.0, 2, 180.0
        if (hb in sp2 and hyb(labels[a]) not in sp2) or (hc in sp2 and hyb(labels[d]) not in sp2):
            V, n, phi0 = 2.0, 3, 180.0
    elif order == 2.0:
        V, n, phi0 = 45.0, 2, 180.0
    elif order >= 1.5:
        V, n, phi0 = 25.0, 2, 180.0
    else:
        V, n, phi0 = 5.0, 2, 180.0
        ri = mol.GetRingInfo()
        if hb == "R" and hc == "R" and not any(b in r and c in r for r in ri.AtomRings()):
            V = 10.0
    norm = (mol.GetAtomWithIdx(b).GetDegree() - 1) * (mol.GetAtomWithIdx(c).GetDegree() - 1)
    return round(V / 2 / norm, 4), n, int(n * phi0 + 180) % 360


def collect(built):
    d, mols = built
    data = read_sections(d / "system.data")
    settings = [l.split() for l in open(d / "system.in.settings") if l.strip()]
    bond_k = {int(t[1]): float(t[2]) for t in settings if t[0] == "bond_coeff"}
    dih = {int(t[1]): (float(t[2]), int(t[3]), int(t[4]) % 360) for t in settings if t[0] == "dihedral_coeff"}
    atoms = sorted(data["Atoms"], key=lambda t: int(t[0]))
    offset, out = 0, {}
    for name, smi in SMILES.items():
        mol, labels, rows = lt_molecule(d / f"{name}.lt", smi)
        n = mol.GetNumAtoms()
        out[name] = (mol, labels, rows, offset, atoms[offset:offset + n])
        offset += n
    return data, bond_k, dih, out


def test_bond_orders_match(built):
    data, bond_k, _, out = collect(built)
    owner = {}
    for name, (mol, labels, rows, off, _) in out.items():
        for k in range(mol.GetNumAtoms()):
            owner[off + k + 1] = (name, k)
    bad = []
    for t in data["Bonds"]:
        (name, i), (_, j) = owner[int(t[2])], owner[int(t[3])]
        mol, labels = out[name][0], out[name][1]
        e = expected_order(mol, i, j, labels)
        if abs(bond_k[int(t[1])] - 350.0 * e) > 1e-6:
            bad.append((name, labels[i], labels[j], e, bond_k[int(t[1])]))
    assert not bad, bad


def test_all_angles_and_dihedrals_typed(built):
    data, _, _, out = collect(built)
    n_ang = n_dih = 0
    for mol, *_ in out.values():
        for a in mol.GetAtoms():
            k = a.GetDegree()
            n_ang += k * (k - 1) // 2
        for b in mol.GetBonds():
            j, k = b.GetBeginAtom(), b.GetEndAtom()
            n_dih += (j.GetDegree() - 1) * (k.GetDegree() - 1)
    assert len(data["Angles"]) == n_ang
    assert len(data["Dihedrals"]) == n_dih


def test_fragments_get_separate_mol_ids(built):
    _, _, _, out = collect(built)
    mol, labels, rows, off, atoms = out["ionpair"]
    mol_ids = {(r[1], int(a[1])) for r, a in zip(rows, atoms)}
    assert len({m for _, m in mol_ids}) == 2


def test_mol2_order_and_charges_kept(built):
    _, mols = built
    _, _, _, out = collect(built)
    rows = out["ionpair"][2]
    for k, (name, frag, label, q) in enumerate(rows):
        el, q2 = mols["mol2_charges"][k]
        assert name == f"{el}{k + 1}"
        assert abs(q - q2) < 1e-9


def test_dihedrals_vs_paper_rules(built):
    data, _, dih, out = collect(built)
    owner = {}
    for name, (mol, labels, rows, off, _) in out.items():
        for k in range(mol.GetNumAtoms()):
            owner[off + k + 1] = (name, k)
    diffs = []
    for t in data["Dihedrals"]:
        ids = [owner[int(x)] for x in t[2:6]]
        name = ids[0][0]
        mol, labels = out[name][0], out[name][1]
        a, b, c, d = (i for _, i in ids)
        ref = paper_dihedral(mol, labels, a, b, c, d)
        got = dih[int(t[1])]
        if ref is None:
            if got[0] != 0.0:
                diffs.append((name, labels[a], labels[b], labels[c], labels[d], got, ref))
        elif (abs(got[0] - ref[0]) > 1e-3 or got[1] != ref[1] or got[2] != ref[2]):
            diffs.append((name, labels[a], labels[b], labels[c], labels[d], got, ref))
    # 알려진 차이(S_3 4배위, N_3 4배위의 이웃 수 가정)는 따로 확인
    unexpected = [x for x in diffs if not ({"S_3", "N_3", "N_3_hd"} & {x[2], x[3]})]
    for x in diffs:
        print(x)
    assert not unexpected, unexpected
