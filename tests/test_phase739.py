"""Phase 739: 9-methylcarbazole — drop indicated-H when N9 is substituted.

IUPAC 2013: when N at the indicated-H position bears a substituent, the
nH prefix is dropped from the parent name.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N9-substituted: no 9H in name
    ("Cn1c2ccccc2c2ccccc21",       "9-methylcarbazole"),
    # C-substituted: 9H stays
    ("Cc1cccc2c1[nH]c1ccccc12",    "1-methyl-9H-carbazole"),
    ("Cc1ccc2c(c1)[nH]c1ccccc12",  "2-methyl-9H-carbazole"),
    ("Cc1ccc2[nH]c3ccccc3c2c1",    "3-methyl-9H-carbazole"),
    ("Cc1cccc2[nH]c3ccccc3c12",    "4-methyl-9H-carbazole"),
    # Regression: unsubstituted parent unchanged
    ("c1ccc2c(c1)[nH]c1ccccc12",   "9H-carbazole"),
])
def test_phase739_carbazole_indicated_h(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
