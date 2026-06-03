"""Phase 344: E/Z in disulfide, sulfoxide, and sulfone parent chains (IUPAC 2013).

Unsaturated parent chains in disulfide, sulfoxide, and sulfone substitutive
names now carry E/Z stereo descriptors and double-bond locants.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # disulfide E/Z
    ("C/C=C/CSSC",              "[(2E)-but-2-en-1-yl] methyl disulfide"),
    (r"C/C=C\CSSC",             "[(2Z)-but-2-en-1-yl] methyl disulfide"),
    # sulfoxide E/Z
    ("C/C=C/CS(=O)C",           "[(2E)-but-2-en-1-yl] methyl sulfoxide"),
    (r"C/C=C\CS(=O)C",          "[(2Z)-but-2-en-1-yl] methyl sulfoxide"),
    # sulfone E/Z
    ("C/C=C/CS(=O)(=O)C",       "[(2E)-but-2-en-1-yl] methyl sulfone"),
    (r"C/C=C\CS(=O)(=O)C",      "[(2Z)-but-2-en-1-yl] methyl sulfone"),
    # regressions: saturated chains unchanged
    ("CSSC",                    "dimethyl disulfide"),
    ("CSSCC",                   "ethyl methyl disulfide"),
    ("CS(=O)C",                 "dimethyl sulfoxide"),
    ("CS(=O)(=O)C",             "dimethyl sulfone"),
])
def test_phase344_ez_disulfide_sulfoxide_sulfone_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
