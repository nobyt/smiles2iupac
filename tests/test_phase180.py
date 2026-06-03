"""Phase 180: アセトアルデヒド・アセトン保留名 (IUPAC 2013 P-31.1.3)

  CC=O    → acetaldehyde  (not ethanal)
  CC(=O)C → acetone       (not propan-2-one)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 保留名
    ("CC=O",         "acetaldehyde"),
    ("CC(=O)C",      "acetone"),
    # 回帰: 他のアルデヒド・ケトンは変わらない
    ("CCC=O",        "propanal"),
    ("CCCC=O",       "butanal"),
    ("CC(=O)CC",     "butan-2-one"),
    ("CCC(=O)CC",    "pentan-3-one"),
    # 回帰: 置換アセトアルデヒドは通常命名 (ethanal ベース)
    ("ClCC=O",       "2-chloroethanal"),
    ("OCC=O",        "2-hydroxyethanal"),
])
def test_phase180_acetaldehyde_acetone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
