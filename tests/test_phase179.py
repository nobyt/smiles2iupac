"""Phase 179: ハロ酢酸・フェニル酢酸保留名 (IUPAC 2013 P-65.1.1.4)

  ClCH2COOH → chloroacetic acid
  BrCH2COOH → bromoacetic acid
  FCH2COOH  → fluoroacetic acid
  ICH2COOH  → iodoacetic acid
  CCl3COOH  → trichloroacetic acid
  PhCH2COOH → phenylacetic acid
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # モノハロ酢酸
    ("OC(=O)CCl",           "chloroacetic acid"),
    ("OC(=O)CBr",           "bromoacetic acid"),
    ("OC(=O)CF",            "fluoroacetic acid"),
    ("OC(=O)CI",            "iodoacetic acid"),
    # ポリハロ酢酸
    ("OC(=O)C(Cl)(Cl)Cl",   "trichloroacetic acid"),
    ("OC(=O)C(Cl)Cl",       "dichloroacetic acid"),
    ("OC(=O)C(F)(F)F",      "trifluoroacetic acid"),
    ("OC(=O)C(F)F",         "difluoroacetic acid"),
    # フェニル酢酸
    ("OC(=O)Cc1ccccc1",     "phenylacetic acid"),
    # 回帰: 置換 > 2C 以上は通常命名
    ("OC(=O)CCCl",          "3-chloropropanoic acid"),
    ("OC(=O)CCC(=O)O",      "succinic acid"),
])
def test_phase179_haloacetic_and_phenylacetic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
