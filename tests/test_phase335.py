"""Phase 335: E/Z in alkenyl substituent names for carbamates, carbamic acids, etc. (IUPAC 2013).

Alkenyl substituents with defined E/Z geometry now include the stereo descriptor
when used as substituent prefixes (e.g., "(2E)-but-2-en-1-yl" in carbamate names).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carbamate with E/Z alkenyl ester group
    ("C/C=C/COC(=O)N",            "(2E)-but-2-en-1-yl carbamate"),
    (r"C/C=C\COC(=O)N",           "(2Z)-but-2-en-1-yl carbamate"),
    # carbamic acid with N-alkenyl substituent
    ("C/C=C/CNC(=O)O",            "N-((2E)-but-2-en-1-yl)carbamic acid"),
    (r"C/C=C\CNC(=O)O",           "N-((2Z)-but-2-en-1-yl)carbamic acid"),
    # regressions: saturated / non-stereo alkenyl unchanged
    ("CC=CCNC(=O)O",              "N-but-2-en-1-ylcarbamic acid"),
    ("CCNC(=O)O",                 "N-ethylcarbamic acid"),
    ("CNC(=O)O",                  "N-methylcarbamic acid"),
    ("CCOC(=O)N",                 "ethyl carbamate"),
    ("COC(=O)N",                  "methyl carbamate"),
])
def test_phase335_ez_alkenyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
