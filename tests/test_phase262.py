"""Phase 262: carboxyoxy substituent — no double parentheses bug (IUPAC 2013).

  OC(=O)CCCOC(=O)O  → 4-(carboxyoxy)butanoic acid
  OC(=O)CCOC(=O)O   → 3-(carboxyoxy)propanoic acid
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carboxyoxy substituent (previously produced double parentheses)
    ("OC(=O)CCCOC(=O)O",  "4-(carboxyoxy)butanoic acid"),
    ("OC(=O)CCOC(=O)O",   "3-(carboxyoxy)propanoic acid"),
    ("OC(=O)COC(=O)O",    "2-(carboxyoxy)acetic acid"),
    # regression: plain carbonate (no chain) unchanged
    ("OC(=O)OCC",         "ethyl hydrogen carbonate"),
])
def test_phase262_carboxyoxy(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
