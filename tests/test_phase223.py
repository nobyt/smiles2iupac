"""Phase 223: selenide, diselenide, telluride, ditelluride (IUPAC 2013 P-63.6).

R-Se-R' → dialkyl selenide; R-Se-Se-R' → dialkyl diselenide.
Tellurium analogs follow the same pattern.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenide
    ("C[Se]C",    "(methylselanyl)methane"),
    ("CC[Se]CC",  "(ethylselanyl)ethane"),
    ("C[Se]CC",   "(methylselanyl)ethane"),
    # diselenide
    ("C[Se][Se]C",    "(methyldiselanyl)methane"),
    ("CC[Se][Se]CC",  "(ethyldiselanyl)ethane"),
    # telluride
    ("C[Te]C",    "(methyltellanyl)methane"),
    # ditelluride
    ("C[Te][Te]C",    "(methylditellanyl)methane"),
    # regression: sulfide and disulfide still work
    ("CSC",       "(methylsulfanyl)methane"),
    ("CSSC",      "(methyldisulfanyl)methane"),
])
def test_phase223_selenide_diselenide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
