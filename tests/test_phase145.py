"""Phase 145: ホスフェート・ホスホネート・ホスフィネートエステル (IUPAC 2013 P-62.4)

trimethyl phosphate, triethyl phosphate,
dimethyl methylphosphonate, diethyl ethylphosphonate,
methyl dimethylphosphinate
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ホスフェートエステル (トリエステル)
    ("COP(=O)(OC)OC",      "trimethyl phosphate"),
    ("CCOP(=O)(OCC)OCC",   "triethyl phosphate"),
    ("CCCOP(=O)(OCCC)OCCC", "tripropyl phosphate"),
    # ホスホネートエステル
    ("COP(=O)(OC)C",       "dimethyl methylphosphonate"),
    ("CCOP(=O)(OCC)CC",    "diethyl ethylphosphonate"),
    ("COP(=O)(OC)CC",      "dimethyl ethylphosphonate"),
    # ホスフィネートエステル
    ("COP(=O)(C)C",        "methyl dimethylphosphinate"),
    # 回帰: 既存リン命名
    ("CP(=O)(O)O",         "methylphosphonic acid"),
    ("CP(C)C",             "trimethylphosphane"),
    ("O=P(O)(O)O",         "phosphoric acid"),
    # 回帰: 一般化合物
    ("CC",                 "ethane"),
    ("CS(=O)C",            "(methylsulfinyl)methane"),
])
def test_phase145_phosphate_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
