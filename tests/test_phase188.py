"""Phase 188: シクロアルキルカルボン酸エステル

  O=C(OC)C1CCCC1  → methyl cyclopentanecarboxylate
  O=C(OC)C1CCCCC1 → methyl cyclohexanecarboxylate
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # シクロアルキルエステル
    ("O=C(OC)C1CCCC1",    "methyl cyclopentanecarboxylate"),
    ("O=C(OC)C1CCCCC1",   "methyl cyclohexanecarboxylate"),
    ("O=C(OC)C1CC1",      "methyl cyclopropanecarboxylate"),
    ("O=C(OC)C1CCC1",     "methyl cyclobutanecarboxylate"),
    ("O=C(OCC)C1CCCCC1",  "ethyl cyclohexanecarboxylate"),
    ("O=C(OCC)C1CCCC1",   "ethyl cyclopentanecarboxylate"),
    # 回帰: 直鎖エステルは変わらない
    ("CC(=O)OC",           "methyl acetate"),
    ("CCCC(=O)OC",         "methyl butanoate"),
    ("CC(=O)OCC",          "ethyl acetate"),
])
def test_phase188_cycloalkyl_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
