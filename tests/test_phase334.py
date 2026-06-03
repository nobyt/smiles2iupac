"""Phase 334: E/Z in sulfonyl halide chains (IUPAC 2013).

Unsaturated sulfonyl halide chains now carry E/Z stereodescriptors.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfonyl chloride E/Z
    ("C/C=C/CS(=O)(=O)Cl",  "(2E)-but-2-ene-1-sulfonyl chloride"),
    (r"C/C=C\CS(=O)(=O)Cl", "(2Z)-but-2-ene-1-sulfonyl chloride"),
    # sulfonyl fluoride E/Z
    ("C/C=C/CS(=O)(=O)F",   "(2E)-but-2-ene-1-sulfonyl fluoride"),
    # regressions: saturated sulfonyl halides unchanged
    ("CCS(=O)(=O)Cl",       "ethanesulfonyl chloride"),
    ("CCCS(=O)(=O)Cl",      "propane-1-sulfonyl chloride"),
    ("c1ccc(S(=O)(=O)Cl)cc1", "benzenesulfonyl chloride"),
])
def test_phase334_sulfonyl_halide_ez(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
