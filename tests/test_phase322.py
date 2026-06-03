"""Phase 322: E/Z in N-substituted amide and thioic acid names (IUPAC 2013 P-93.5).

N-substituted amides and thioic acids of alpha,beta-unsaturated acids
now carry the E/Z prefix. Also fixes thioic acid suffix ('enthioic' → 'enethioic').
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-substituted amide + E/Z
    ("C/C=C/C(=O)NC",               "(2E)-N-methylbut-2-enamide"),
    ("C/C=C/C(=O)NCC",              "(2E)-N-ethylbut-2-enamide"),
    (r"C/C=C\C(=O)NC",              "(2Z)-N-methylbut-2-enamide"),
    ("C/C=C/C(=O)N(C)C",            "(2E)-N,N-dimethylbut-2-enamide"),
    # thioic S-acid + E/Z
    ("C/C=C/C(=O)S",                "(2E)-but-2-enethioic S-acid"),
    (r"C/C=C\C(=O)S",               "(2Z)-but-2-enethioic S-acid"),
    # thioic O-acid + E/Z
    ("C/C=C/C(=S)O",                "(2E)-but-2-enethioic O-acid"),
    # thioic acid enethioic suffix fix (no stereo markers)
    ("CC=CC(=O)S",                  "but-2-enethioic S-acid"),
    # regressions: saturated N-sub amides unchanged
    ("CC(=O)NC",                    "N-methylacetamide"),
    ("CC(=O)N(C)C",                 "N,N-dimethylacetamide"),
    # regressions: saturated thioic acids unchanged
    ("CC(=O)S",                     "ethanethioic S-acid"),
    ("CC(=S)O",                     "ethanethioic O-acid"),
])
def test_phase322_ez_amide_thioic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
