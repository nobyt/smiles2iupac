"""Phase 681: acridin-9(10H)-one and quinoxalin-2(1H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acridin-9(10H)-one: C(9,=O)-N(10,H)-benzo-benzo (C2v-symmetric: 1≡8,2≡7,3≡6,4≡5)
    # (parent in Phase 413; C9=O not methylable; C1-C4 methylable via C2v)
    ("O=c1c2ccccc2[nH]c2ccccc12",    "acridin-9(10H)-one"),
    ("O=c1c2c(C)cccc2[nH]c2ccccc12", "1-methylacridin-9(10H)-one"),
    ("O=c1c2cc(C)ccc2[nH]c2ccccc12", "2-methylacridin-9(10H)-one"),
    ("O=c1c2ccc(C)cc2[nH]c2ccccc12", "3-methylacridin-9(10H)-one"),
    ("O=c1c2cccc(C)c2[nH]c2ccccc12", "4-methylacridin-9(10H)-one"),
    # quinoxalin-2(1H)-one: N(1,H)-C(2,=O)-C(3)-N(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 413; C2=O and N4 not methylable; C3 and C5-C8 methylable)
    ("O=c1cnc2ccccc2[nH]1",    "quinoxalin-2(1H)-one"),
    ("O=c1c(C)nc2ccccc2[nH]1", "3-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2c(C)cccc2[nH]1", "5-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2cc(C)ccc2[nH]1", "6-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2ccc(C)cc2[nH]1", "7-methylquinoxalin-2(1H)-one"),
    ("O=c1cnc2cccc(C)c2[nH]1", "8-methylquinoxalin-2(1H)-one"),
])
def test_phase681(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
