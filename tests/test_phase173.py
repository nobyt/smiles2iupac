"""Phase 173: 環付き酸ハライド命名 (IUPAC 2013 P-65.1.2.1)

脂環式環に結合したカルボニルハライドを cycloalkanecarbonyl halide 形式で命名:
  O=C(Cl)C1CCCCC1  → cyclohexanecarbonyl chloride
  O=C(Cl)C1CCCC1   → cyclopentanecarbonyl chloride
  O=C(Br)C1CCCCC1  → cyclohexanecarbonyl bromide
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # シクロアルカン結合酸ハライド
    ("O=C(Cl)C1CCCCC1",  "cyclohexanecarbonyl chloride"),
    ("O=C(Cl)C1CCCC1",   "cyclopentanecarbonyl chloride"),
    ("O=C(Cl)C1CCC1",    "cyclobutanecarbonyl chloride"),
    ("O=C(Cl)C1CC1",     "cyclopropanecarbonyl chloride"),
    ("O=C(Br)C1CCCCC1",  "cyclohexanecarbonyl bromide"),
    ("O=C(F)C1CCCCC1",   "cyclohexanecarbonyl fluoride"),
    # 置換環
    ("O=C(Cl)C1CCCCC1C", "2-methylcyclohexanecarbonyl chloride"),
    # 回帰: 鎖状酸ハライドは変わらない
    ("CC(=O)Cl",          "acetyl chloride"),
    ("CCCC(=O)Cl",        "butanoyl chloride"),
    # 芳香環 (benzoyl) は変わらない
    ("O=C(Cl)c1ccccc1",   "benzoyl chloride"),
])
def test_phase173_cyclic_acid_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
