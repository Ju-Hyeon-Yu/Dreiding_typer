# CLI 입출력(조각 배치, mol2 대응 오류)을 moltemplate 없이 확인하는 테스트
import os
import sys

import numpy as np
import pytest
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dreiding_lt  # noqa: E402


def test_ion_fragments_do_not_overlap():
    m = Chem.AddHs(Chem.MolFromSmiles("Cn1cc[n+](Cc2ccc(I)cc2)c1.[O-]S(=O)(=O)C(F)(F)F.[I-]"))
    xyz = np.array(dreiding_lt.embed(m, seed=1))
    frags = Chem.GetMolFrags(m)
    for a in range(len(frags)):
        for b in range(a + 1, len(frags)):
            d = np.linalg.norm(xyz[list(frags[a])][:, None] - xyz[list(frags[b])][None], axis=2)
            assert d.min() > 3.0


def test_mol2_mismatch_exits(tmp_path):
    mol2 = tmp_path / "x.mol2"
    mol2.write_text("@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 MOL 0.0\n2 O2 1.4 0 0 O.3 1 MOL 0.0\n"
                    "@<TRIPOS>BOND\n1 1 2 1\n")
    with pytest.raises(SystemExit):
        dreiding_lt.main(["--smiles", "CCO", "--mol2", str(mol2), "-o", str(tmp_path / "x.lt")])


@pytest.mark.parametrize("name, atype, el", [
    ("C1", "C.ar", "C"), ("CL3", "Cl", "Cl"), ("C7", "ca", "C"), ("Cl2", "cl", "Cl"),
    ("I1", "i", "I"), ("H12", "hc", "H"), ("N4", "N.pl3", "N"), ("X1", "Br", "Br"),
])
def test_mol2_element_detection(name, atype, el):
    assert dreiding_lt._element(name, atype) == el
