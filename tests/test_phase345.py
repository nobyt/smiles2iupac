"""Phase 345: E/Z in selenide, telluride, diselenide, ditelluride parent chains (IUPAC 2013).

Unsaturated parent chains in selenide, telluride, and their di-chalcogen
analogues now carry E/Z descriptors and double-bond locants in the
substitutive name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenide E/Z
    ("C/C=C/C[Se]C",            "(2E)-1-(methylselanyl)but-2-ene"),
    (r"C/C=C\C[Se]C",           "(2Z)-1-(methylselanyl)but-2-ene"),
    # telluride E/Z
    ("C/C=C/C[Te]C",            "(2E)-1-(methyltellanyl)but-2-ene"),
    (r"C/C=C\C[Te]C",           "(2Z)-1-(methyltellanyl)but-2-ene"),
    # regressions: saturated chains unchanged
    ("C[Se]C",                  "(methylselanyl)methane"),
    ("C[Se]CC",                 "(methylselanyl)ethane"),
    ("C[Te]C",                  "(methyltellanyl)methane"),
])
def test_phase345_ez_selenide_telluride_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
