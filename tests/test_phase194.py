"""Phase 194: アシルアジド命名 (IUPAC 2013 P-66.6)

  CC(=O)N=[N+]=[N-]   → acetyl azide
  CCC(=O)N=[N+]=[N-]  → propanoyl azide

構造: C(=O)-N=[N+]=[N-] — アジド基を持つアシル化合物。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # アシルアジド
    ("C(=O)N=[N+]=[N-]",       "formyl azide"),
    ("CC(=O)N=[N+]=[N-]",      "acetyl azide"),
    ("CCC(=O)N=[N+]=[N-]",     "propanoyl azide"),
    ("CCCC(=O)N=[N+]=[N-]",    "butanoyl azide"),
    # 回帰: アミドは変わらない
    ("CC(=O)N",                "acetamide"),
    ("CC(=O)NC",               "N-methylacetamide"),
    ("CC(=O)NN",               "ethanohydrazide"),
    # 回帰: アルキルアジド (PIN: 置換命名)
    ("[N-]=[N+]=NC",           "azidomethane"),
    ("[N-]=[N+]=NCC",          "azidoethane"),
    # 回帰: エステルは変わらない
    ("CC(=O)OC",               "methyl acetate"),
])
def test_phase194_acyl_azide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
