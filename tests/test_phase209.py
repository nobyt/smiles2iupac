"""Phase 209: acyl peroxide naming (IUPAC 2013 P-65.1.5.3)

  CC(=O)OOC(=O)C    → diethanoyl peroxide
  CCC(=O)OOC(=O)CC  → dipropanoyl peroxide
  CC(=O)OOC(=O)CC   → ethanoyl propanoyl peroxide (mixed)

RC(=O)-O-O-C(=O)R' are named di(acyl) peroxides (symmetric) or
acyl1 acyl2 peroxide (asymmetric).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric acyl peroxides
    ("CC(=O)OOC(=O)C",    "diethanoyl peroxide"),
    ("CCC(=O)OOC(=O)CC",  "dipropanoyl peroxide"),
    ("C(=O)OOC=O",        "diformyl peroxide"),
    # Asymmetric (alphabetical order)
    ("CC(=O)OOC(=O)CC",   "ethanoyl propanoyl peroxide"),
    # regression: simple peroxide unaffected
    ("COOC",              "dimethyl peroxide"),
    ("CCOOCC",            "diethyl peroxide"),
    # regression: peroxyacid unaffected
    ("CC(=O)OO",          "ethaneperoxoic acid"),
])
def test_phase209_acyl_peroxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
