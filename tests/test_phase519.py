"""Phase 519: selenoxide / selenone naming (IUPAC 2013 P-65.3.2)

C[Se](=O)C → dimethyl selenoxide; C[Se](=O)(=O)C → dimethyl selenone.
Previously both were mis-detected as selenide.

Phase 856: the tellurium analogs (C[Te](=O)C, C[Te](=O)(=O)C) had the same
bug -- mis-detected as plain "telluride" (a copy-paste error: functional_group.py
had `"selenone" if is_se else "telluride"` instead of `else "tellurone"`,
and similarly for selenoxide/telluroxide), silently dropping the oxo/dioxo
information entirely.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenoxide: C-Se(=O)-C
    ("C[Se](=O)C",       "dimethyl selenoxide"),
    ("CC[Se](=O)CC",     "diethyl selenoxide"),
    ("CC[Se](=O)C",      "ethyl methyl selenoxide"),
    # selenone: C-Se(=O)₂-C
    ("C[Se](=O)(=O)C",   "dimethyl selenone"),
    ("CC[Se](=O)(=O)CC", "diethyl selenone"),
    # telluroxide: C-Te(=O)-C (Phase 856)
    ("C[Te](=O)C",       "dimethyl telluroxide"),
    ("CC[Te](=O)CC",     "diethyl telluroxide"),
    ("CC[Te](=O)C",      "ethyl methyl telluroxide"),
    # tellurone: C-Te(=O)₂-C (Phase 856)
    ("C[Te](=O)(=O)C",   "dimethyl tellurone"),
    ("CC[Te](=O)(=O)CC", "diethyl tellurone"),
    # regression: plain selenide/telluride unaffected
    ("C[Se]C",           "methylselanylmethane"),
    ("CC[SeH]",          "ethaneselenol"),
    ("C[Te]C",           "methyltellanylmethane"),
    # regression: sulfoxide/sulfone unaffected
    ("CS(=O)C",          "dimethyl sulfoxide"),
    ("CS(=O)(=O)C",      "dimethyl sulfone"),
])
def test_phase519_selenoxide_selenone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
