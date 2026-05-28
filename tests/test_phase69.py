"""Phase 69: シアン酸エステル (O-C≡N) / チオシアン酸エステル (S-C≡N)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # シアン酸エステル (cyanate)
    ("COC#N", "methyl cyanate"),
    ("CCOC#N", "ethyl cyanate"),
    ("CCCOC#N", "propyl cyanate"),
    # チオシアン酸エステル (thiocyanate)
    ("CSC#N", "methyl thiocyanate"),
    ("CCSC#N", "ethyl thiocyanate"),
    ("CCCSC#N", "propyl thiocyanate"),
])
def test_phase69_cyanate_thiocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
