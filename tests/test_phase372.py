"""Phase 372: Phosphonothioate ester naming (R-P(=S)(OR')2).

Previously R-P(=S)(OR')2 fell through all phosphorus detection branches
(which required P=O) and produced the wrong name '(P)methane'. Added
detection for the P=S analog of phosphonate esters and a new naming
function that uses 'phosphonothioate' as the suffix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric ester groups
    ("CP(=S)(OC)OC",    "dimethyl methylphosphonothioate"),
    ("CP(=S)(OCC)OCC",  "diethyl methylphosphonothioate"),
    # Asymmetric ester groups would require more testing
    # Different P-alkyl
    ("CCP(=S)(OC)OC",   "dimethyl ethylphosphonothioate"),
    # Regressions: phosphonate ester (P=O) unchanged
    ("CP(=O)(OC)OC",    "dimethyl methylphosphonate"),
    ("CCP(=O)(OCC)OCC", "diethyl ethylphosphonate"),
    # Other P compounds unchanged
    ("CB(O)O",          "methylboronic acid"),
])
def test_phase372_phosphonothioate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
