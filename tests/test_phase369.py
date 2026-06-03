"""Phase 369: Sulfone/sulfoxide preferred IUPAC 2013 names (dialkyl sulfone format).

Previously _name_sulfoxide_sulfone used the substitutive form
'(methylsulfonyl)methane' / '(methylsulfinyl)methane'. IUPAC 2013 P-65.3.1
requires the functional-class form: group names + 'sulfone'/'sulfoxide'.

Rules:
  - Same groups: di<group> sulfone/sulfoxide (e.g. dimethyl sulfone)
  - Different groups: <group1> <group2> sulfone/sulfoxide (alphabetical)
  - Complex group names (containing digits) → wrap in parentheses
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric simple
    ("CS(=O)(=O)C",           "dimethyl sulfone"),
    ("CS(=O)C",               "dimethyl sulfoxide"),
    ("CCS(=O)(=O)CC",         "diethyl sulfone"),
    ("CCS(=O)CC",             "diethyl sulfoxide"),
    # Asymmetric simple (alphabetical)
    ("CS(=O)(=O)CC",          "ethyl methyl sulfone"),
    ("CS(=O)CC",              "ethyl methyl sulfoxide"),
    # Aryl groups
    ("c1ccccc1S(=O)(=O)c1ccccc1",  "diphenyl sulfone"),
    ("c1ccccc1S(=O)c1ccccc1",      "diphenyl sulfoxide"),
    # Vinyl groups
    ("C=CS(=O)(=O)C",         "ethenyl methyl sulfone"),
    ("C=CS(=O)C",             "ethenyl methyl sulfoxide"),
    # Complex group: prop-1-en-1-yl needs parens
    ("CC=CS(=O)(=O)C",        "methyl (prop-1-en-1-yl) sulfone"),
    # Complex with E/Z: but-2-en-1-yl needs brackets
    ("C/C=C/CS(=O)(=O)C",     "[(2E)-but-2-en-1-yl] methyl sulfone"),
    # Regressions: sulfonamide, sulfonic acid, sulfonate not affected
    ("CS(=O)(=O)N",           "methanesulfonamide"),
    ("CS(=O)(=O)O",           "methanesulfonic acid"),
    ("CS(=O)(=O)OC",          "methyl methanesulfonate"),
])
def test_phase369_sulfone_sulfoxide_preferred(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
