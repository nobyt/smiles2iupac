"""Phase 268: benzophenone retained name; benzene dicarboxylate esters; methyl salicylate (IUPAC 2013).

  O=C(c1ccccc1)c1ccccc1       → benzophenone               (retained PIN)
  CCOC(=O)c1ccccc1C(=O)OCC   → diethyl benzene-1,2-dicarboxylate
  COC(=O)c1ccccc1O            → methyl salicylate           (retained)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzophenone (retained PIN, IUPAC 2013 P-31.1.3.4)
    ("O=C(c1ccccc1)c1ccccc1",            "benzophenone"),
    # benzene dicarboxylic acid diesters (phthalate, isophthalate, terephthalate)
    ("CCOC(=O)c1ccccc1C(=O)OCC",         "diethyl benzene-1,2-dicarboxylate"),
    ("CCOC(=O)c1cccc(C(=O)OCC)c1",       "diethyl benzene-1,3-dicarboxylate"),
    ("CCOC(=O)c1ccc(C(=O)OCC)cc1",       "diethyl benzene-1,4-dicarboxylate"),
    ("COC(=O)c1ccccc1C(=O)OC",           "dimethyl benzene-1,2-dicarboxylate"),
    # methyl salicylate (methyl ester of salicylic acid, IUPAC 2013 retained)
    ("COC(=O)c1ccccc1O",                  "methyl salicylate"),
    # regression: acetophenone and simple ketones unchanged
    ("O=C(C)c1ccccc1",                    "acetophenone"),
    ("O=C(CC)c1ccccc1",                   "1-phenylpropan-1-one"),
    # regression: mixed aryl ketones unchanged
    ("O=C(c1ccccc1)c1ccc(Cl)cc1",        "(4-chlorophenyl)(phenyl)methanone"),
    # regression: simple diester unchanged
    ("CCOC(=O)CC(=O)OCC",                "diethyl malonate"),
])
def test_phase268_benzophenone_dicarboxylate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
