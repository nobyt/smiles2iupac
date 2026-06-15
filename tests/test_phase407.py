"""Phase 407: Isoindolin-1-one and N-substituted lactam retained names.

IUPAC 2013 P-31.1.3: isoindolin-1-one is the retained name for the
5-membered lactam ring fused with benzene (O at C1, N at position 2).

Also fixes 1-substituted indolin-2-one (oxindole) naming by ensuring
the exo C=O is part of the base retained name rather than an oxo prefix.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoindolin-1-one (N-H)
    ("O=C1NCc2ccccc21",               "isoindolin-1-one"),
    # N-methylisoindolin-1-one
    ("O=C1N(C)Cc2ccccc21",            "2-methylisoindolin-1-one"),
    # 1-methylindolin-2-one
    ("CN1C(=O)Cc2ccccc21",            "1-methylindolin-2-one"),
    # regression: indolin-2-one unchanged
    ("O=C1Cc2ccccc2N1",               "indolin-2-one"),
    # regression: 1-methylindoline unchanged
    ("CN1CCc2ccccc21",                "1-methylindoline"),
    # regression: pyrrolidine-2,5-dione unchanged
    ("O=C1CCC(=O)N1",                 "pyrrolidine-2,5-dione"),
    # regression: phthalimide unchanged
    ("O=C1NC(=O)c2ccccc21",           "phthalimide"),
    # regression: benzene unchanged
    ("c1ccccc1",                       "benzene"),
])
def test_phase407_isoindolin_1_one(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
