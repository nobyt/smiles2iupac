"""Phase 147: 複数成分 SMILES (塩命名) & ジカルボキシレート (IUPAC 2013)

salt naming from dot-notation SMILES:
  sodium acetate, disodium sulfate, dipotassium succinate,
  azanium acetate, sodium bicarbonate

dicarboxylate anion naming:
  oxalate, malonate, succinate, glutarate, adipate
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジカルボキシレートアニオン
    ("[O-]C(=O)C(=O)[O-]",             "oxalate"),
    ("[O-]C(=O)CC(=O)[O-]",             "malonate"),
    ("[O-]C(=O)CCC(=O)[O-]",            "succinate"),
    ("[O-]C(=O)CCCC(=O)[O-]",           "glutarate"),
    ("[O-]C(=O)CCCCC(=O)[O-]",          "adipate"),
    # モノカルボキシレート
    ("O=C[O-]",                          "formate"),
    ("CC(=O)[O-]",                       "acetate"),
    ("CCC(=O)[O-]",                      "propanoate"),
    # 無機炭酸系アニオン
    ("O=C([O-])O",                       "bicarbonate"),
    ("O=C([O-])[O-]",                    "carbonate"),
    # 塩命名 (dot notation)
    ("[Na+].[O-]C(=O)C",                 "sodium acetate"),
    ("[Na+].[Na+].[O-]S(=O)(=O)[O-]",   "disodium sulfate"),
    ("[K+].[K+].[O-]C(=O)CCC(=O)[O-]",  "dipotassium succinate"),
    ("[NH4+].CC(=O)[O-]",                "azanium acetate"),
    ("[Na+].[O-]C(=O)O",                 "sodium bicarbonate"),
    ("[Ca+2].[O-]C(=O)C.[O-]C(=O)C",    "calcium diacetate"),
    # 回帰: 通常カルボン酸
    ("CC(=O)O",                          "acetic acid"),
    ("OC(=O)CCC(=O)O",                   "succinic acid"),
    ("CC",                               "ethane"),
])
def test_phase147_salts_and_dianions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
