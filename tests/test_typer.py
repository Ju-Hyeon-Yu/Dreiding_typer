# 대표 분자별 DREIDING 라벨 할당 결과를 확인하는 단위 테스트
import os
import sys

import pytest
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dreiding_typer import DreidingFF, TypingError, assign_labels  # noqa: E402

FF = DreidingFF()


def typed(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    labs, warns = assign_labels(mol, FF)
    heavy = [labs[a.GetIdx()] for a in mol.GetAtoms() if a.GetSymbol() != "H"]
    hyd = sorted({labs[a.GetIdx()] for a in mol.GetAtoms() if a.GetSymbol() == "H"})
    return heavy, hyd, warns, mol, labs


def check_orders(mol, labs):
    """라벨 쌍의 결합 차수가 실제와 같은지 (sp2 계열끼리만).

    두 원자가 모두 R 이면 1.5, 방향족 결합인데 R 이 아니면(비방향족으로 바꾼 피리돈류 고리) 케쿨레 차수,
    나머지는 SMILES 차수.
    """
    kek = Chem.Mol(mol)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        sp2 = all(any(t in labs[k] for t in ("_2", "_R", "_1")) for k in (i, j))
        if not sp2:
            continue
        if "_R" in labs[i] and "_R" in labs[j] and b.GetIsAromatic():
            e = 1.5
        else:
            e = kek.GetBondWithIdx(b.GetIdx()).GetBondTypeAsDouble()
        assert FF.bond_order(labs[i], labs[j]) == e, (labs[i], labs[j], e)


@pytest.mark.parametrize("smiles, heavy, hyd", [
    ("c1ccccc1", ["C_R"] * 6, ["H"]),
    ("c1ccncc1", ["C_R", "C_R", "C_R", "N_R_d1_ha", "C_R", "C_R"], ["H"]),
    ("c1cnc[nH]1", ["C_R", "C_R", "N_R_d1_ha", "C_R", "N_R_d2_hd"], ["H", "H_HB"]),
    ("Cn1cc[n+](C)c1", ["C_3", "N_R_d2", "C_R", "C_R", "N_R_d2", "C_3", "C_R"], ["H"]),
    ("c1ccoc1", ["C_R", "C_R", "C_R", "O_R_ha", "C_R"], ["H"]),          # moltemplate 푸란 예제 (+_ha)
    ("CC(N)=O", ["C_3", "C_2_b1", "N_2_b1_d2_hd", "O_2_ha"], ["H", "H_HB"]),
    ("CC(=O)[O-]", ["C_3", "C_2_b1", "O_2_ha", "O_2_b1_ha"], ["H"]),
    ("CC(=O)O", ["C_3", "C_2", "O_2_ha", "O_3_hd"], ["H", "H_HB"]),
    ("Nc1ccccc1", ["N_2_d2_hd"] + ["C_R"] * 6, ["H", "H_HB"]),
    ("NC(N)=[NH2+]", ["N_2_d2_hd", "C_2_b2", "N_2_d2_hd", "N_2_b2_d2_hd"], ["H_HB"]),
    ("CS(N)(=O)=O", ["C_3", "S_3", "N_3_hd", "O_2_ha", "O_2_ha"], ["H", "H_HB"]),
    ("F[B-](F)(F)F", ["F_ha", "B_3", "F_ha", "F_ha", "F_ha"], []),
    ("Oc1ccccc1", ["O_3_hd"] + ["C_R"] * 6, ["H", "H_HB"]),
    ("Ic1ccccc1", ["I"] + ["C_R"] * 6, ["H"]),
    ("C[N+](C)(C)C", ["C_3", "N_3", "C_3", "C_3", "C_3"], ["H"]),
    ("C[NH3+]", ["C_3", "N_3_hd"], ["H", "H_HB"]),
    ("CC(C)=N", ["C_3", "C_2", "C_3", "N_2_d1_hd"], ["H", "H_HB"]),     # 매뉴얼 그림 2c 와 같음
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
    assert h[:4] == ["O_2_ha", "N_2_b1_d2", "O_2_b1_ha", "C_R"]
    check_orders(mol, labs)


def test_amide_follows_manual_text_rule():
    """매뉴얼 그림 2b 는 N 에만 _b1 이 있어 C–N 이 이중결합이 되므로, 본문 규칙(짝)대로 C 에도 _b1."""
    h, _, _, mol, labs = typed("CC(N)=O")
    assert h[1:3] == ["C_2_b1", "N_2_b1_d2_hd"]
    assert FF.bond_order("C_2_b1", "N_2_b1_d2_hd") == 1.0
    assert FF.bond_order("C_2", "N_2_b1_d2_hd") == 2.0      # 그림 그대로일 때
    check_orders(mol, labs)


def test_conjugated_chains_satisfy_orders():
    for smi in ["C=CC=CC=C", "O=C1C=CC=C1", "C=CC(=O)O", "O=C(O)c1ccccc1", "CC(=O)Nc1ccccc1",
                "[O-]C(=O)c1ccccc1", "NC(N)=O", "CN(C)C=C", "c1ccn(c1)-c1ccccc1"]:
        _, _, _, mol, labs = typed(smi)
        check_orders(mol, labs)


SMILES_SET = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "web", "tests", "smiles_set.txt")


def test_manual_rules_hold_on_smiles_set():
    """R 은 방향족 고리 원자에만(황은 S_3 뿐이라 예외), _b1/_b2 는 같은 플래그 이웃과 짝, 결합 차수 일치."""
    n = 0
    for smi in open(SMILES_SET).read().split():
        try:
            _, _, _, mol, labs = typed(smi)
        except TypingError:
            continue
        n += 1
        kek = Chem.Mol(mol)
        Chem.Kekulize(kek, clearAromaticFlags=True)
        for a in mol.GetAtoms():
            lab = labs[a.GetIdx()]
            if "_R" in lab:
                assert a.GetIsAromatic(), (smi, a.GetIdx(), lab)
            for fl in ("_b1", "_b2"):
                if fl in lab:
                    assert any(fl in labs[x.GetIdx()] for x in a.GetNeighbors()), (smi, a.GetIdx(), lab)
        check_orders(mol, labs)
    assert n > 100


@pytest.mark.parametrize("smiles", ["O=c1cccc[nH]1", "O=c1cc[nH]c(=O)[nH]1", "O=c1ccc2ccccc2o1",
                                    "Cn1cnc2c1c(=O)n(C)c(=O)n2C"])
def test_exocyclic_carbonyl_ring_becomes_nonaromatic(smiles):
    h, _, warns, mol, labs = typed(smiles)
    for b in mol.GetBonds():
        if b.GetBondTypeAsDouble() == 2.0 and not b.GetIsAromatic():
            c = b.GetBeginAtom() if b.GetBeginAtom().GetSymbol() == "C" else b.GetEndAtom()
            assert labs[c.GetIdx()].startswith("C_2"), (smiles, labs[c.GetIdx()])
    assert any("사이클로펜타디에논" in w for w in warns)
    check_orders(mol, labs)


def test_coumarin_keeps_benzene_ring_aromatic():
    h, _, _, _, _ = typed("O=c1ccc2ccccc2o1")
    assert h.count("C_R") == 6 and h[-1] == "O_3_ha"


def test_triflate_warns_on_sulfur():
    h, _, warns, _, _ = typed("[O-]S(=O)(=O)C(F)(F)F")
    assert h == ["O_2_ha", "S_3", "O_2_ha", "O_2_ha", "C_3", "F_ha", "F_ha", "F_ha"]
    assert any("S" in w and "S_3" in w for w in warns)


def test_tfsi_imide_nitrogen():
    h, _, warns, _, _ = typed("FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F")
    assert "N_3_ha" in h
    assert any("2배위 N" in w for w in warns)


def test_no_force_field_value_warning():
    _, _, warns, _, _ = typed("CCCl")
    assert not any("질량" in w for w in warns)


@pytest.mark.parametrize("smiles", ["[O-]Cl(=O)(=O)=O", "c1ccc([I+]c2ccccc2)cc1"])
def test_hypervalent_halogen_raises(smiles):
    with pytest.raises(TypingError, match="초원자가"):
        typed(smiles)


def test_unsupported_element_raises():
    with pytest.raises(TypingError):
        typed("[Se]=C=[Se]")
