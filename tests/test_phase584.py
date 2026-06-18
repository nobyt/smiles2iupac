"""Phase 584: Substituted 1H-benzo[f]indole and 1H-benzo[g]indole naming.
N-H at position 1; substitutable C positions: 2-9 (8 positions each).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-benzo[f]indole
    ("c1ccc2cc3[nH]ccc3cc2c1",       "1H-benzo[f]indole"),
    ("Cc1cc2cc3ccccc3cc2[nH]1",      "2-methyl-1H-benzo[f]indole"),
    ("Cc1c[nH]c2cc3ccccc3cc12",      "3-methyl-1H-benzo[f]indole"),
    ("Cc1c2ccccc2cc2[nH]ccc12",      "4-methyl-1H-benzo[f]indole"),
    ("Cc1cccc2cc3[nH]ccc3cc12",      "5-methyl-1H-benzo[f]indole"),
    ("Cc1ccc2cc3[nH]ccc3cc2c1",      "6-methyl-1H-benzo[f]indole"),
    ("Cc1ccc2cc3cc[nH]c3cc2c1",      "7-methyl-1H-benzo[f]indole"),
    ("Cc1cccc2cc3cc[nH]c3cc12",      "8-methyl-1H-benzo[f]indole"),
    ("Cc1c2ccccc2cc2cc[nH]c12",      "9-methyl-1H-benzo[f]indole"),
    # 1H-benzo[g]indole
    ("c1ccc2c(c1)ccc1cc[nH]c12",     "1H-benzo[g]indole"),
    ("Cc1cc2ccc3ccccc3c2[nH]1",      "2-methyl-1H-benzo[g]indole"),
    ("Cc1c[nH]c2c1ccc1ccccc12",      "3-methyl-1H-benzo[g]indole"),
    ("Cc1cc2ccccc2c2[nH]ccc12",      "4-methyl-1H-benzo[g]indole"),
    ("Cc1cc2cc[nH]c2c2ccccc12",      "5-methyl-1H-benzo[g]indole"),
    ("Cc1cccc2c1ccc1cc[nH]c12",      "6-methyl-1H-benzo[g]indole"),
    ("Cc1ccc2c(ccc3cc[nH]c32)c1",    "7-methyl-1H-benzo[g]indole"),
    ("Cc1ccc2ccc3cc[nH]c3c2c1",      "8-methyl-1H-benzo[g]indole"),
    ("Cc1cccc2ccc3cc[nH]c3c12",      "9-methyl-1H-benzo[g]indole"),
])
def test_phase584_benzo_fg_indole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
