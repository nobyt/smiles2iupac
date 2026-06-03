"""Phase 321: E/Z stereodescriptor in acid anhydride names (IUPAC 2013 P-93.5).

Anhydrides of alpha,beta-unsaturated acids now carry the E/Z prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # symmetric E anhydride
    ("C/C=C/C(=O)OC(=O)/C=C/C",    "(2E)-but-2-enoic anhydride"),
    # symmetric Z anhydride
    (r"C/C=C\C(=O)OC(=O)/C=C\C",   "(2Z)-but-2-enoic anhydride"),
    # mixed E/Z anhydride
    (r"C/C=C/C(=O)OC(=O)/C=C\C",   "(2E)-but-2-enoic (2Z)-but-2-enoic anhydride"),
    # regressions: saturated anhydrides unchanged
    ("CC(=O)OC(=O)C",               "acetic anhydride"),
    ("CCCC(=O)OC(=O)CCC",          "butanoic anhydride"),
    ("CC(=O)OC(=O)CCC",            "acetic butanoic anhydride"),
])
def test_phase321_ez_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
