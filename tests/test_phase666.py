"""Phase 666: quinoxalin-2(1H)-one and quinazoline-2,4(1H,3H)-dione methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinoxalin-2(1H)-one: N(1,H)-C(2,=O)-C(3)-N(4)-C4a-...-C8a
    # (parent in Phase 413; C2=O, N1-H, aromatic N4 not methylable)
    ("O=c1cnc2ccccc2[nH]1",    "quinoxalin-2(1H)-one"),
    ("O=c1c(C)nc2ccccc2[nH]1", "3-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2c(C)cccc2[nH]1", "5-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2cc(C)ccc2[nH]1", "6-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2ccc(C)cc2[nH]1", "7-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2cccc(C)c2[nH]1", "8-methylquinoxalin-2(1H)-one"),
    # quinazoline-2,4(1H,3H)-dione: N(1,H)-C(2,=O)-N(3,H)-C(4,=O)-C4a-...-C8a
    # (parent in Phase 415; C2=O, C4=O, N1-H, N3-H not methylable; C2-symmetric: 5≡8, 6≡7)
    ("O=c1[nH]c(=O)c2ccccc2[nH]1",    "quinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2c(C)cccc2[nH]1", "5-methylquinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2cc(C)ccc2[nH]1", "6-methylquinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2ccc(C)cc2[nH]1", "7-methylquinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2cccc(C)c2[nH]1", "8-methylquinazoline-2,4(1H,3H)-dione"),
])
def test_phase666(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
