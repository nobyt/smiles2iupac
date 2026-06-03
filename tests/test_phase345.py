"""Phase 345: E/Z in selenide, telluride, diselenide, ditelluride parent chains (IUPAC 2013).

Unsaturated parent chains in selenide, telluride, and their di-chalcogen
analogues now carry E/Z descriptors and double-bond locants in the
substitutive name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenide E/Z
    ("C/C=C/C[Se]C",            "[(2E)-but-2-en-1-yl] methyl selenide"),
    (r"C/C=C\C[Se]C",           "[(2Z)-but-2-en-1-yl] methyl selenide"),
    # telluride E/Z
    ("C/C=C/C[Te]C",            "[(2E)-but-2-en-1-yl] methyl telluride"),
    (r"C/C=C\C[Te]C",           "[(2Z)-but-2-en-1-yl] methyl telluride"),
    # regressions: saturated chains unchanged
    ("C[Se]C",                  "dimethyl selenide"),
    ("C[Se]CC",                 "ethyl methyl selenide"),
    ("C[Te]C",                  "dimethyl telluride"),
])
def test_phase345_ez_selenide_telluride_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
