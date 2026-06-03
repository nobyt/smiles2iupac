"""Phase 190: アシルオキシ置換基命名 (acetyloxy, butanoyloxy など)

  CC(=O)OCC(=O)O   → 2-(acetyloxy)acetic acid
  CCCC(=O)OCCC(=O)O → 3-(butanoyloxy)propanoic acid

IUPAC P-65.1.2.3: アシルオキシ基は {acyl}oxy として命名する。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # カルボン酸主官能基 + アシルオキシ置換基
    ("CC(=O)OCC(=O)O",       "2-(acetyloxy)acetic acid"),
    ("CCCC(=O)OCCC(=O)O",    "3-(butanoyloxy)propanoic acid"),
    ("CCC(=O)OCC(=O)O",      "2-(propanoyloxy)acetic acid"),
    ("CC(=O)OCCC(=O)O",      "3-(acetyloxy)propanoic acid"),
    # 回帰: 通常エーテル酸素 (methoxy/ethoxy) は変わらない
    ("COCC(=O)O",             "2-methoxyacetic acid"),
    ("CCOCC(=O)O",            "2-ethoxyacetic acid"),
    ("COC",                   "methoxymethane"),
    ("CCOC",                  "methoxyethane"),
    # 回帰: エステル命名は変わらない
    ("CC(=O)OC",              "methyl acetate"),
    ("CCCC(=O)OCC",           "ethyl butanoate"),
])
def test_phase190_acyloxy_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
