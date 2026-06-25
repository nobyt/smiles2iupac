"""Phase 667: acridin-9(10H)-one methyl derivatives.

C2v-symmetric (C1≡C8, C2≡C7, C3≡C6, C4≡C5); lowest-locant preferred.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acridin-9(10H)-one: C1-C4 ring A / C9(=O) / N10(H) / C5-C8 ring B
    # (parent in Phase 413; C9=O and N10-H not methylable; 4 unique positions)
    ("O=c1c2ccccc2[nH]c2ccccc12",    "acridin-9(10H)-one"),
    ("O=c1c2c(C)cccc2[nH]c2ccccc12", "1-methylacridin-9(10H)-one"),
    ("O=c1c2cc(C)ccc2[nH]c2ccccc12", "2-methylacridin-9(10H)-one"),
    ("O=c1c2ccc(C)cc2[nH]c2ccccc12", "3-methylacridin-9(10H)-one"),
    ("O=c1c2cccc(C)c2[nH]c2ccccc12", "4-methylacridin-9(10H)-one"),
])
def test_phase667(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
