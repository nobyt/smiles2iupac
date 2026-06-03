"""Phase 344: E/Z in disulfide, sulfoxide, and sulfone parent chains (IUPAC 2013).

Unsaturated parent chains in disulfide, sulfoxide, and sulfone substitutive
names now carry E/Z stereo descriptors and double-bond locants.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # disulfide E/Z
    ("C/C=C/CSSC",              "(2E)-1-(methyldisulfanyl)but-2-ene"),
    (r"C/C=C\CSSC",             "(2Z)-1-(methyldisulfanyl)but-2-ene"),
    # sulfoxide E/Z
    ("C/C=C/CS(=O)C",           "(2E)-1-(methylsulfinyl)but-2-ene"),
    (r"C/C=C\CS(=O)C",          "(2Z)-1-(methylsulfinyl)but-2-ene"),
    # sulfone E/Z
    ("C/C=C/CS(=O)(=O)C",       "(2E)-1-(methylsulfonyl)but-2-ene"),
    (r"C/C=C\CS(=O)(=O)C",      "(2Z)-1-(methylsulfonyl)but-2-ene"),
    # regressions: saturated chains unchanged
    ("CSSC",                    "(methyldisulfanyl)methane"),
    ("CSSCC",                   "(methyldisulfanyl)ethane"),
    ("CS(=O)C",                 "(methylsulfinyl)methane"),
    ("CS(=O)(=O)C",             "(methylsulfonyl)methane"),
])
def test_phase344_ez_disulfide_sulfoxide_sulfone_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
