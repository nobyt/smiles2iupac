"""Phase 180: アセトアルデヒド保留名; アセトン → PIN: propan-2-one (IUPAC 2013 P-65.1.2.2)

  CC=O    → acetaldehyde  (not ethanal; retained PIN)
  CC(=O)C → propan-2-one  (acetone は保留名; PIN は propan-2-one)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acetaldehyde は保留 PIN; acetone は保留名 → PIN へ
    ("CC=O",         "acetaldehyde"),
    ("CC(=O)C",      "propan-2-one"),
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
