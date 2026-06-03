"""Phase 332: E/Z in sulfonic, sulfinic, sulfenic, sulfonamide, and sulfinamide (IUPAC 2013).

Unsaturated aliphatic sulfonic/sulfinic/sulfenic acids and sulfonamide/sulfinamide
chains now carry E/Z stereodescriptors.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfonic acid E/Z
    ("C/C=C/CS(=O)(=O)O",    "(2E)-but-2-ene-1-sulfonic acid"),
    (r"C/C=C\CS(=O)(=O)O",   "(2Z)-but-2-ene-1-sulfonic acid"),
    # sulfinic acid E/Z
    ("C/C=C/CS(=O)O",        "(2E)-but-2-ene-1-sulfinic acid"),
    (r"C/C=C\CS(=O)O",       "(2Z)-but-2-ene-1-sulfinic acid"),
    # sulfenic acid E/Z
    ("C/C=C/CSO",            "(2E)-but-2-ene-1-sulfenic acid"),
    # sulfonamide E/Z
    ("C/C=C/CS(=O)(=O)N",    "(2E)-but-2-ene-1-sulfonamide"),
    (r"C/C=C\CS(=O)(=O)N",   "(2Z)-but-2-ene-1-sulfonamide"),
    # sulfinamide E/Z
    ("C/C=C/CS(=O)N",        "(2E)-but-2-ene-1-sulfinamide"),
    (r"C/C=C\CS(=O)N",       "(2Z)-but-2-ene-1-sulfinamide"),
    # regressions: saturated chains unchanged
    ("CCS(=O)(=O)O",         "ethanesulfonic acid"),
    ("CCCS(=O)O",            "propane-1-sulfinic acid"),
    ("CS(=O)(=O)N",          "methanesulfonamide"),
    ("CS(=O)(=O)NC",         "N-methylmethanesulfonamide"),
])
def test_phase332_sulfur_acid_amide_ez(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
