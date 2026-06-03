"""Phase 338: E/Z parenthesization in N-alkenyl substituents across all N-handlers (IUPAC 2013).

When an alkenyl N-substituent carries an E/Z stereo descriptor such as
"(2E)-but-2-en-1-yl", IUPAC requires an additional outer parenthesis layer
to avoid ambiguity: N-[(2E)-but-2-en-1-yl]…, not N-(2E)-but-2-en-1-yl….
This phase fixes that for: imine, nitrone-imine, hydroxylamine, selenoamide,
telluramide, sulfamic acid, sulfinamide, carbamate, and aryl amine handlers.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imine N-substituent
    ("C/C=C/CN=C",           "N-[(2E)-but-2-en-1-yl]methanimine"),
    (r"C/C=C\CN=C",          "N-[(2Z)-but-2-en-1-yl]methanimine"),
    # hydroxylamine N-substituent
    ("C/C=C/CNO",            "N-[(2E)-but-2-en-1-yl]hydroxylamine"),
    (r"C/C=C\CNO",           "N-[(2Z)-but-2-en-1-yl]hydroxylamine"),
    # selenoamide N-substituent
    ("C/C=C/CNC(=[Se])C",    "N-[(2E)-but-2-en-1-yl]ethaneselenoamide"),
    # telluramide N-substituent
    ("C/C=C/CNC(=[Te])C",    "N-[(2E)-but-2-en-1-yl]ethaneteluramide"),
    # sulfamic acid N-substituent
    ("C/C=C/CNS(=O)(=O)O",   "N-[(2E)-but-2-en-1-yl]sulfamic acid"),
    # sulfinamide N-substituent
    ("C/C=C/CNS(=O)C",       "N-[(2E)-but-2-en-1-yl]methanesulfinamide"),
    # carbamate N-substituent
    ("C/C=C/CNC(=O)OC",      "methyl N-[(2E)-but-2-en-1-yl]carbamate"),
    # aniline N-substituent
    ("C/C=C/CNc1ccccc1",     "N-[(2E)-but-2-en-1-yl]aniline"),
    (r"C/C=C\CNc1ccccc1",    "N-[(2Z)-but-2-en-1-yl]aniline"),
    # regressions: saturated N-substituents unchanged
    ("CNC=C",                "N-methylethenamine"),
    ("CNO",                  "N-methylhydroxylamine"),
    ("CCNO",                 "N-ethylhydroxylamine"),
    ("CNC(=[Se])C",          "N-methylethaneselenoamide"),
])
def test_phase338_ez_n_alkenyl_parenthesization(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
