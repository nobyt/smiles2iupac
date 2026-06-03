"""Phase 337: E/Z in secondary/tertiary amine parent chains (IUPAC 2013).

When the principal chain of a secondary or tertiary amine contains a C=C double
bond with defined E/Z geometry, the stereo descriptor must appear in the name.
N-substituents that are alkenyl with E/Z descriptors also need parenthesization.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # E/Z on the parent chain of secondary amine
    ("C/C=C/CNC",               "(2E)-N-methylbut-2-en-1-amine"),
    (r"C/C=C\CNC",              "(2Z)-N-methylbut-2-en-1-amine"),
    ("CNC/C=C/C",               "(2E)-N-methylbut-2-en-1-amine"),
    # E/Z parent chain with ethyl N-substituent (4C alkenyl beats 2C ethyl)
    ("C/C=C/CNCC",              "(2E)-N-ethylbut-2-en-1-amine"),
    (r"C/C=C\CNCC",             "(2Z)-N-ethylbut-2-en-1-amine"),
    # E/Z parent chain tertiary amine
    ("C/C=C/CN(C)C",            "(2E)-N,N-dimethylbut-2-en-1-amine"),
    # E/Z on N-substituent (alkenyl is shorter than saturated parent)
    ("CCCCCNC/C=C/C",           "N-[(2E)-but-2-en-1-yl]pentan-1-amine"),
    (r"CCCCCNC/C=C\C",          "N-[(2Z)-but-2-en-1-yl]pentan-1-amine"),
    # regressions: saturated secondary/tertiary amines unchanged
    ("CNC",                     "N-methylmethanamine"),
    ("CCNCC",                   "N-ethylethanamine"),
    ("CN(C)C",                  "N,N-dimethylmethanamine"),
    ("CCNC",                    "N-methylethanamine"),
])
def test_phase337_ez_amine_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
