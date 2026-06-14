"""Phase 521: cyanate/thiocyanate functional-class PIN (IUPAC 2013 P-65.3.1)

R-O-C≡N → {R} cyanate; R-S-C≡N → {R} thiocyanate.
Same functional-class pattern as isocyanate/isothiocyanate (Phase 516).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # cyanate
    ("COC#N",              "methyl cyanate"),
    ("CCOC#N",             "ethyl cyanate"),
    ("CCCOC#N",            "propyl cyanate"),
    ("c1ccccc1OC#N",       "phenyl cyanate"),
    # thiocyanate
    ("CSC#N",              "methyl thiocyanate"),
    ("CCSC#N",             "ethyl thiocyanate"),
    ("CCCSC#N",            "propyl thiocyanate"),
    ("N#CSc1ccccc1",       "phenyl thiocyanate"),
    # acids regression
    ("OC#N",               "cyanic acid"),
    ("SC#N",               "thiocyanic acid"),
    # isocyanate/isothiocyanate regression
    ("CN=C=O",             "methyl isocyanate"),
    ("CN=C=S",             "methyl isothiocyanate"),
])
def test_phase521_cyanate_thiocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
