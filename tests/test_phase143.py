"""Phase 143: リン・ホウ素・ケイ素化合物命名 (IUPAC 2013 P-62, P-67, P-68)

phosphane, diphosphane, phosphoric acid, phosphorous acid,
borane, silane, disilane,
alkylphosphanes, phosphonic acid, boronic acid, borinic acid,
organic boranes, organic silanes
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 無機リン化合物 (保留名)
    ("P",            "phosphane"),
    ("PP",           "diphosphane"),
    ("O=P(O)(O)O",   "phosphoric acid"),
    ("OP(O)O",       "phosphorous acid"),
    ("O=[PH](O)O",   "hypophosphorous acid"),
    # アルキルホスファン
    ("CP",           "methylphosphane"),
    ("CPC",          "dimethylphosphane"),
    ("CP(C)C",       "trimethylphosphane"),
    ("CCP",          "ethylphosphane"),
    # ホスホン酸・ホスフィン酸
    ("CP(=O)(O)O",   "methylphosphonic acid"),
    ("CCP(=O)(O)O",  "ethylphosphonic acid"),
    ("CP(=O)(O)C",   "dimethylphosphinic acid"),
    # 無機ホウ素化合物 (保留名)
    ("B",            "borane"),
    # ボロン酸・ボリン酸
    ("CB(O)O",       "methylboronic acid"),
    ("CCB(O)O",      "ethylboronic acid"),
    ("OB(O)c1ccccc1", "phenylboronic acid"),
    ("CB(O)C",       "dimethylborinic acid"),
    # 有機ボラン
    ("CB(C)C",       "trimethylborane"),
    # 無機ケイ素化合物 (保留名)
    ("[SiH4]",       "silane"),
    ("[SiH3][SiH3]", "disilane"),
    # 有機シラン
    ("C[SiH3]",      "methylsilane"),
    ("C[SiH2]C",     "dimethylsilane"),
    ("C[SiH](C)C",   "trimethylsilane"),
    ("C[Si](C)(C)C", "tetramethylsilane"),
    # 回帰: 通常の有機化合物
    ("CC",           "ethane"),
    ("CCCC",         "butane"),
    ("CS(=O)C",      "(methylsulfinyl)methane"),
    ("c1ccncc1",     "pyridine"),
])
def test_phase143_pbs_compounds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
