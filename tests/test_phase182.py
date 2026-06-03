"""Phase 182: ケト酸保留名 (IUPAC 2013 P-65.1.1.4)

  CC(=O)C(=O)O   → pyruvic acid    (not 2-oxopropanoic acid)
  CC(=O)CCC(=O)O → levulinic acid  (not 4-oxopentanoic acid)
  CC(=O)CC(=O)O  → acetoacetic acid (not 3-oxobutanoic acid)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 保留名
    ("CC(=O)C(=O)O",    "pyruvic acid"),
    ("CC(=O)CCC(=O)O",  "levulinic acid"),
    ("CC(=O)CC(=O)O",   "acetoacetic acid"),
    # 回帰: 他のケト酸は通常命名
    ("CCC(=O)C(=O)O",   "2-oxobutanoic acid"),
    ("O=CCC(=O)C",      "3-oxobutanal"),
    ("O=CCCC(=O)C",     "4-oxopentanal"),
])
def test_phase182_keto_acid_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
