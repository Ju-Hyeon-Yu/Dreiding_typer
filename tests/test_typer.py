# 대표 분자별 DREIDING 라벨 할당 결과를 확인하는 단위 테스트
import os
import sys

import pytest
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dreiding_typer import TypingError, assign_labels, implied_order, load_labels  # noqa: E402

LABELS = load_labels()


def typed(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    labs, warns = assign_labels(mol, LABELS)
    heavy = [labs[a.GetIdx()] for a in mol.GetAtoms() if a.GetSymbol() != "H"]
    hyd = sorted({labs[a.GetIdx()] for a in mol.GetAtoms() if a.GetSymbol() == "H"})
    return heavy, hyd, warns, mol, labs


def check_orders(mol, labs):
    """방향족 1.5, 나머지 결합 차수와 라벨 차수가 같은지 (sp2 계열끼리만)."""
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        sp2 = all(any(t in labs[k] for t in ("_2", "_R", "_1")) for k in (i, j))
        if not sp2:
            continue
        resonant = b.GetIsAromatic() or ("_R" in labs[i] and "_R" in labs[j] and
                                         not (b.GetBeginAtom().IsInRing() or b.GetEndAtom().IsInRing()))
        e = 1.5 if resonant else b.GetBondTypeAsDouble()
        assert implied_order(labs[i], labs[j]) == e, (labs[i], labs[j], e)


@pytest.mark.parametrize("smiles, heavy, hyd", [
    ("c1ccccc1", ["C_R"] * 6, ["H"]),
    ("c1ccncc1", ["C_R", "C_R", "C_R", "N_R_d1_ha", "C_R", "C_R"], ["H"]),
    ("c1cnc[nH]1", ["C_R", "C_R", "N_R_d1_ha", "C_R", "N_R_d2_hd"], ["H", "H_HB"]),
    ("Cn1cc[n+](C)c1", ["C_3", "N_R_d2", "C_R", "C_R", "N_R_d2", "C_3", "C_R"], ["H"]),
    ("CC(N)=O", ["C_3", "C_2", "N_R_d2_hd", "O_2_ha"], ["H", "H_HB"]),
    ("CC(=O)[O-]", ["C_3", "C_R", "O_R_ha", "O_R_ha"], ["H"]),
    ("CC(=O)O", ["C_3", "C_2", "O_2_ha", "O_3_hd"], ["H", "H_HB"]),
    ("Oc1ccccc1", ["O_3_hd"] + ["C_R"] * 6, ["H", "H_HB"]),
    ("Ic1ccccc1", ["I"] + ["C_R"] * 6, ["H"]),
    ("C[N+](C)(C)C", ["C_3", "N_3", "C_3", "C_3", "C_3"], ["H"]),
    ("C[NH3+]", ["C_3", "N_3_hd"], ["H", "H_HB"]),
    ("CC(C)=N", ["C_3", "C_2", "C_3", "N_2_d1_hd"], ["H", "H_HB"]),
    ("CC#N", ["C_3", "C_1", "N_1_ha"], ["H"]),
    ("CCO", ["C_3", "C_3", "O_3_hd"], ["H", "H_HB"]),
    ("COC", ["C_3", "O_3_ha", "C_3"], ["H"]),
])
def test_labels(smiles, heavy, hyd):
    h, hy, _, mol, labs = typed(smiles)
    assert h == heavy
    assert hy == hyd
    check_orders(mol, labs)


def test_butadiene_matches_manual():
    # 라벨 매뉴얼 그림 1a
    h, _, _, _, _ = typed("C=CC=C")
    assert h == ["C_2", "C_2_b1", "C_2_b1", "C_2"]


def test_biphenyl_bridge_b1():
    h, _, _, mol, labs = typed("c1ccc(cc1)-c1ccccc1")
    assert h.count("C_R_b1") == 2
    check_orders(mol, labs)


def test_nitrobenzene():
    h, _, _, mol, labs = typed("O=[N+]([O-])c1ccccc1")
    assert h[:4] == ["O_R_ha", "N_R_b1_d2", "O_R_ha", "C_R_b1"]
    check_orders(mol, labs)


def test_conjugated_chains_satisfy_orders():
    for smi in ["C=CC=CC=C", "O=C1C=CC=C1", "C=CC(=O)O", "O=C(O)c1ccccc1", "CC(=O)Nc1ccccc1"]:
        _, _, _, mol, labs = typed(smi)
        check_orders(mol, labs)


def test_triflate_warns_on_sulfur():
    h, _, warns, _, _ = typed("[O-]S(=O)(=O)C(F)(F)F")
    assert h == ["O_2_ha", "S_3", "O_2_ha", "O_2_ha", "C_3", "F_ha", "F_ha", "F_ha"]
    assert any("S" in w and "S_3" in w for w in warns)


def test_tfsi_imide_nitrogen():
    h, _, _, _, _ = typed("FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F")
    assert "N_R_d1_ha" in h


def test_unsupported_element_raises():
    with pytest.raises(TypingError):
        typed("[Se]=C=[Se]")
