"""Phase 223: selenide, diselenide, telluride, ditelluride (IUPAC 2013 P-63.6).

R-Se-R' → dialkyl selenide; R-Se-Se-R' → dialkyl diselenide.
Tellurium analogs follow the same pattern.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenide
    ("C[Se]C",    "dimethyl selenide"),
    ("CC[Se]CC",  "diethyl selenide"),
    ("C[Se]CC",   "ethyl methyl selenide"),
    # diselenide
    ("C[Se][Se]C",    "dimethyl diselenide"),
    ("CC[Se][Se]CC",  "diethyl diselenide"),
    # telluride
    ("C[Te]C",    "dimethyl telluride"),
    # ditelluride
    ("C[Te][Te]C",    "dimethyl ditelluride"),
    # regression: sulfide and disulfide still work
    ("CSC",       "dimethyl sulfide"),
    ("CSSC",      "dimethyl disulfide"),
])
def test_phase223_selenide_diselenide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
