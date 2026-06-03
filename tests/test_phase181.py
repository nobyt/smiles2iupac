"""Phase 181: アクリル酸系保留名 (IUPAC 2013 P-65.1.1.4)

  C=CC(=O)O    → acrylic acid      (not prop-2-enoic acid)
  CC(=C)C(=O)O → methacrylic acid  (not 2-methylprop-2-enoic acid)
  OCC(=O)O     → glycolic acid     (not 2-hydroxyacetic acid)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 保留名
    ("C=CC(=O)O",       "acrylic acid"),
    ("CC(=C)C(=O)O",    "methacrylic acid"),
    ("OCC(=O)O",        "glycolic acid"),
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
