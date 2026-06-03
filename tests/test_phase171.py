"""Phase 171: フェニル置換基内の括弧修正 + Phase 170 テスト統合

フェニル上の置換基名がロカントを含む場合に括弧が付くよう修正:
  CC(C)Cc1ccc(cc1)C(C)C(=O)O  → 2-(4-(2-methylpropyl)phenyl)propanoic acid
  CC(C)Cc1ccc(N)cc1            → 4-(2-methylpropyl)aniline
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ロカント付きフェニル置換基への括弧
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O",  "2-[4-(2-methylpropyl)phenyl]propanoic acid"),
    ("CC(C)Cc1ccc(N)cc1",             "4-(2-methylpropyl)aniline"),
    ("CC(C)Cc1ccc(O)cc1",             "4-(2-methylpropyl)phenol"),
    # ロカントなし置換基: 括弧なし
    ("Cc1ccc(Cl)cc1",                 "1-chloro-4-methylbenzene"),
    ("CCc1ccccc1",                    "ethylbenzene"),
    # Phase 170: 単純無機化合物 retained names
    ("OO",    "hydrogen peroxide"),
    ("[H][H]", "dihydrogen"),
    ("S",      "hydrogen sulfide"),
    ("O",      "water"),
])
def test_phase171_phenyl_brackets_and_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
