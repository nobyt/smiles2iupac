"""Phase 677: 2H-pyran-2-one and 4H-pyran-4-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2H-pyran-2-one (alpha-pyrone): O(1)-C(2,=O)-C(3)-C(4)-C(5)-C(6)
    # (parent in Phase 420; C2=O and O1 not methylable; C3-C6 methylable)
    ("O=c1cccco1",     "2H-pyran-2-one"),
    ("O=c1c(C)ccco1",  "3-methyl-2H-pyran-2-one"),
    ("O=c1cc(C)cco1",  "4-methyl-2H-pyran-2-one"),
    ("O=c1ccc(C)co1",  "5-methyl-2H-pyran-2-one"),
    ("O=c1cccc(C)o1",  "6-methyl-2H-pyran-2-one"),
    # 4H-pyran-4-one (gamma-pyrone): C(2)-C(3)-C(4,=O)-C(5)-C(6)-O(1)
    # (parent in Phase 420; C2v-symmetric: 2≡6, 3≡5; C4=O and O1 not methylable)
    ("O=c1ccocc1",     "4H-pyran-4-one"),
    ("O=c1cc(C)occ1",  "2-methyl-4H-pyran-4-one"),
    ("O=c1c(C)cocc1",  "3-methyl-4H-pyran-4-one"),
])
def test_phase677(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
