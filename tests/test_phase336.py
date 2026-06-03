"""Phase 336: E/Z in N-alkenyl substituents in urea, thiourea, thioamide, and sulfonamide (IUPAC 2013).

N-alkenyl substituents with E/Z stereo now get correct parenthesization in
urea, thiourea, thioamide, and sulfonamide naming.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-alkenyl urea
    ("C/C=C/CNC(=O)N",             "N-[(2E)-but-2-en-1-yl]urea"),
    (r"C/C=C\CNC(=O)N",            "N-[(2Z)-but-2-en-1-yl]urea"),
    # N-alkenyl thiourea
    ("C/C=C/CNC(=S)N",             "N-[(2E)-but-2-en-1-yl]thiourea"),
    # N-alkenyl thioamide
    ("C/C=C/CNC(=S)C",             "N-[(2E)-but-2-en-1-yl]ethanethioamide"),
    # N-alkenyl sulfonamide
    ("C/C=C/CNS(=O)(=O)C",         "N-[(2E)-but-2-en-1-yl]methanesulfonamide"),
    (r"C/C=C\CNS(=O)(=O)C",        "N-[(2Z)-but-2-en-1-yl]methanesulfonamide"),
    # regressions: saturated N-substituents unchanged
    ("CNC(=O)N",                   "N-methylurea"),
    ("CNC(=S)N",                   "N-methylthiourea"),
    ("CS(=O)(=O)NC",               "N-methylmethanesulfonamide"),
    ("CC(=S)NC",                   "N-methylethanethioamide"),
])
def test_phase336_ez_n_alkenyl_urea_thioamide_sulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
