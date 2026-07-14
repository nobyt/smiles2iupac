"""Phase 223: selenide, diselenide, telluride, ditelluride (IUPAC 2013 P-63.6).

R-Se-R' -> substitutive (alkylselanyl)parent (Phase 856); R-Se-Se-R' ->
dialkyl diselenide (functional-class, no substitutive equivalent commonly
used for -Se-Se-). Tellurium analogs follow the same pattern.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenide
    ("C[Se]C",    "methylselanylmethane"),
    ("CC[Se]CC",  "ethylselanylethane"),
    ("C[Se]CC",   "methylselanylethane"),
    # diselenide
    ("C[Se][Se]C",    "dimethyl diselenide"),
    ("CC[Se][Se]CC",  "diethyl diselenide"),
    # telluride
    ("C[Te]C",    "methyltellanylmethane"),
    # ditelluride
    ("C[Te][Te]C",    "dimethyl ditelluride"),
    # regression: sulfide and disulfide still work
    ("CSC",       "(methylsulfanyl)methane"),
    ("CSSC",      "dimethyl disulfide"),
])
def test_phase223_selenide_diselenide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
