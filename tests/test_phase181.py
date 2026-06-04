"""Phase 181: 保留名 (IUPAC 2013 P-65.1.1.4)

  OCC(=O)O     → glycolic acid (retained PIN)

Phase 384 corrected: acrylic acid / methacrylic acid are not IUPAC 2013 PIN:
  C=CC(=O)O    → prop-2-enoic acid
  CC(=C)C(=O)O → 2-methylprop-2-enoic acid
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 保留名 (glycolic acid は PIN)
    ("OCC(=O)O",        "glycolic acid"),
    # acrylic / methacrylic acid: systematic names (Phase 384)
    ("C=CC(=O)O",       "prop-2-enoic acid"),
    ("CC(=C)C(=O)O",    "2-methylprop-2-enoic acid"),
    # 回帰: 置換誘導体は通常命名
    ("ClCC(=O)O",       "chloroacetic acid"),
    ("CC=CC(=O)O",      "but-2-enoic acid"),
    ("OC(=O)CC(=O)O",   "malonic acid"),
    # 回帰: エステルは通常命名
    ("C=CC(=O)OC",      "methyl prop-2-enoate"),
    ("C=CC(=O)OCC",     "ethyl prop-2-enoate"),
])
def test_phase181_acrylic_acid_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
