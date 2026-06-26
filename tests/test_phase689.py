"""Phase 689: dibenz[a,h]anthracene methyl derivatives (positions 1–7; C2 sym so 8–14 duplicate 1–7)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2c(c1)ccc1cc3c(ccc4ccccc43)cc12",  "dibenz[a,h]anthracene"),
    ("Cc1cc2ccccc2c2cc3ccc4ccccc4c3cc12",     "1-methyldibenz[a,h]anthracene"),
    ("Cc1c2ccc3ccccc3c2cc2ccc3ccccc3c12",     "2-methyldibenz[a,h]anthracene"),
    ("Cc1cccc2ccc3cc4c(ccc5ccccc54)cc3c12",   "3-methyldibenz[a,h]anthracene"),
    ("Cc1ccc2ccc3cc4c(ccc5ccccc54)cc3c2c1",   "4-methyldibenz[a,h]anthracene"),
    ("Cc1ccc2c(ccc3cc4c(ccc5ccccc54)cc32)c1", "5-methyldibenz[a,h]anthracene"),
    ("Cc1cccc2c1ccc1cc3c(ccc4ccccc43)cc12",   "6-methyldibenz[a,h]anthracene"),
    ("Cc1cc2cc3c(ccc4ccccc43)cc2c2ccccc12",   "7-methyldibenz[a,h]anthracene"),
])
def test_phase689(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
