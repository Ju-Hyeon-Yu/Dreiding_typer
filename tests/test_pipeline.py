# --lammps 파이프라인(moltemplate.sh + cleanup + 입력 보정)을 LAMMPS 로 읽어 확인하는 테스트
import os
import re
import shutil
import sys

import pytest
from rdkit import Chem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import dreiding_lt  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("moltemplate.sh") is None, reason="moltemplate 없음")

IONPAIR = "Cn1cc[n+](Cc2ccc(I)cc2)c1.[O-]S(=O)(=O)C(F)(F)F"   # 수소결합 주개 없음
PROTIC = "CC(=O)[O-].C[NH3+]"                                 # 주개 있음


def write_mol2(smiles, path):
    """조각을 떨어뜨린 좌표와 형식 전하를 고르게 나눈 전하로 mol2 작성."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    xyz = dreiding_lt.embed(m, seed=3)
    frags = Chem.GetMolFrags(m)
    q = [0.0] * m.GetNumAtoms()
    for fr in frags:
        tot = sum(m.GetAtomWithIdx(i).GetFormalCharge() for i in fr)
        for i in fr:
            q[i] = tot / len(fr)
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\nt\n{m.GetNumAtoms()} {m.GetNumBonds()}\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for a in m.GetAtoms():
            x, y, z = xyz[a.GetIdx()]
            f.write(f"{a.GetIdx() + 1} {a.GetSymbol()}{a.GetIdx() + 1} {x:.4f} {y:.4f} {z:.4f} "
                    f"{a.GetSymbol()} 1 MOL {q[a.GetIdx()]:.8f}\n")
        f.write("@<TRIPOS>BOND\n")
        for k, b in enumerate(m.GetBonds(), 1):
            f.write(f"{k} {b.GetBeginAtomIdx() + 1} {b.GetEndAtomIdx() + 1} 1\n")
    return q


def single_point(d):
    lammps = pytest.importorskip("lammps")
    cwd = os.getcwd()
    os.chdir(d)
    try:
        lmp = lammps.lammps(cmdargs=["-screen", "none", "-log", "none"])
        lmp.file("system.in.init")
        lmp.command("read_data system.data")
        lmp.file("system.in.settings")
        lmp.command("thermo_style custom pe evdwl ecoul elong ebond eangle edihed eimp")
        lmp.command("run 0")
        return {k: lmp.get_thermo(k) for k in ("pe", "evdwl", "ecoul", "elong", "ebond", "eangle", "edihed")}
    finally:
        os.chdir(cwd)


def test_ionpair_periodic_without_donors(tmp_path):
    q = write_mol2(IONPAIR, tmp_path / "ip.mol2")
    run = tmp_path / "run"
    dreiding_lt.main(["--smiles", IONPAIR, "--mol2", str(tmp_path / "ip.mol2"),
                      "-o", str(tmp_path / "ip.lt"), "--lammps", str(run)])
    init = (run / "system.in.init").read_text()
    assert "hbond/dreiding" not in init and "kspace_style" in init
    data = (run / "system.data").read_text()
    assert int(re.search(r"(\d+)\s+atom types", data).group(1)) < 15        # cleanup 적용
    charges = [float(l.split()[3]) for l in re.search(r"\nAtoms[^\n]*\n\n(.*?)\n\n", data, re.S)
               .group(1).splitlines()]
    assert charges == pytest.approx(q, abs=1e-6)
    e = single_point(run)
    assert all(abs(v) < 1e3 for v in e.values()), e


def test_protic_vacuum_with_hbond(tmp_path):
    write_mol2(PROTIC, tmp_path / "p.mol2")
    run = tmp_path / "run"
    dreiding_lt.main(["--smiles", PROTIC, "--mol2", str(tmp_path / "p.mol2"),
                      "-o", str(tmp_path / "p.lt"), "--lammps", str(run), "--vacuum"])
    init = (run / "system.in.init").read_text()
    assert init.startswith("boundary        s s s")
    assert "kspace_style" not in init
    assert re.search(r"hybrid/overlay hbond/dreiding/lj .* lj/cut/coul/cut 100", init)
    assert "lj/cut/coul/long" not in (run / "system.in.settings").read_text()
    e = single_point(run)
    assert e["elong"] == 0.0
    assert all(abs(v) < 1e3 for v in e.values()), e
