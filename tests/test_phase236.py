"""Phase 236: cyclic imide systematic names (IUPAC 2013 P-66.8.3 PIN).

succinimide → pyrrolidine-2,5-dione  (5-membered, IUPAC 2013 PIN)
glutarimide → piperidine-2,6-dione   (6-membered, IUPAC 2013 PIN)
(retained names 'succinimide'/'glutarimide' are not PINs per IUPAC 2013)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=C1CCC(=O)N1",   "pyrrolidine-2,5-dione"),
    ("O=C1NC(=O)CC1",   "pyrrolidine-2,5-dione"),   # same compound, different SMILES
    ("O=C1CCCC(=O)N1",  "piperidine-2,6-dione"),
    ("O=C1NC(=O)CCC1",  "piperidine-2,6-dione"),
    ("O=C1NC(=O)c2ccccc21", "phthalimide"),
    # regression: larger ring uses systematic name
    ("O=C1NC(=O)CCCC1", "azepane-2,7-dione"),
    # regression: lactam (only one C=O) unchanged
    ("O=C1CCCN1",       "pyrrolidin-2-one"),
])
def test_phase236_cyclic_imide_systematic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
