"""Phase 237: additional inorganic retained names (IUPAC 2013).

Selenium oxoacids, pyrophosphoric acid, silicon halides.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=[Se](O)O",           "selenious acid"),
    ("O=[Se](=O)(O)O",       "selenic acid"),
    ("OP(=O)(O)OP(=O)(O)O",  "pyrophosphoric acid"),
    ("[Si](Cl)(Cl)(Cl)Cl",   "tetrachlorosilane"),
    ("[SiH2](Cl)Cl",         "dichlorosilane"),
    ("[SiH3]Cl",             "chlorosilane"),
    # regression: phosphoric acid unchanged
    ("O=P(O)(O)O",           "phosphoric acid"),
    ("O=[Se](O)O",           "selenious acid"),
])
def test_phase237_inorganic_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
