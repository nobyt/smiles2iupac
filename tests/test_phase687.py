"""Phase 687: benz[a]anthracene methyl derivatives (positions 1–12)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benz[a]anthracene: ring D (benz, upper-right) = pos 1-4; ring C = 5,6 (K-region);
    # ring B = 7,12; ring A = 8-11; bay region between pos 1 and 12
    ("c1ccc2cc3c(ccc4ccccc43)cc2c1",   "benz[a]anthracene"),
    ("Cc1cccc2ccc3cc4ccccc4cc3c12",    "1-methylbenz[a]anthracene"),
    ("Cc1ccc2ccc3cc4ccccc4cc3c2c1",    "2-methylbenz[a]anthracene"),
    ("Cc1ccc2c(ccc3cc4ccccc4cc32)c1",  "3-methylbenz[a]anthracene"),
    ("Cc1cccc2c1ccc1cc3ccccc3cc12",    "4-methylbenz[a]anthracene"),
    ("Cc1cc2cc3ccccc3cc2c2ccccc12",    "5-methylbenz[a]anthracene"),
    ("Cc1cc2ccccc2c2cc3ccccc3cc12",    "6-methylbenz[a]anthracene"),
    ("Cc1c2ccccc2cc2c1ccc1ccccc12",    "7-methylbenz[a]anthracene"),
    ("Cc1cccc2cc3c(ccc4ccccc43)cc12",  "8-methylbenz[a]anthracene"),
    ("Cc1ccc2cc3c(ccc4ccccc43)cc2c1",  "9-methylbenz[a]anthracene"),
    ("Cc1ccc2cc3ccc4ccccc4c3cc2c1",    "10-methylbenz[a]anthracene"),
    ("Cc1cccc2cc3ccc4ccccc4c3cc12",    "11-methylbenz[a]anthracene"),
    ("Cc1c2ccccc2cc2ccc3ccccc3c12",    "12-methylbenz[a]anthracene"),
])
def test_phase687(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
