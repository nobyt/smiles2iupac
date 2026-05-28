"""Phase 120: IUPAC 2013 保留名 — formic acid, formamide, acetamide, etc. (P-65.1.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # formic acid (methanoic acid の PIN)
    ("OC=O", "formic acid"),
    ("C(=O)O", "formic acid"),
    # formamide (methanamide の PIN)
    ("NC=O", "formamide"),
    # acetamide (ethanamide の PIN)
    ("CC(=O)N", "acetamide"),
    # formyl chloride (methanoyl chloride の PIN)
    ("ClC=O", "formyl chloride"),
    # acetyl chloride (ethanoyl chloride の PIN)
    ("CC(=O)Cl", "acetyl chloride"),
    # methyl formate / ethyl formate (methanoate の PIN)
    ("COC=O", "methyl formate"),
    ("CCOC=O", "ethyl formate"),
    # acetophenone (1-phenylethan-1-one の PIN)
    ("CC(=O)c1ccccc1", "acetophenone"),
    ("O=C(C)c1ccccc1", "acetophenone"),
    # N-置換体にも保留名を継承
    ("CNC=O", "N-methylformamide"),
    ("CC(=O)NC", "N-methylacetamide"),
    ("CC(=O)N(C)C", "N,N-dimethylacetamide"),
    # 回帰: 3炭素以上は系統名
    ("CCCO", "propan-1-ol"),
    ("CCC(=O)O", "propanoic acid"),
    ("CCC(=O)N", "propanamide"),
    ("CCC(=O)Cl", "propanoyl chloride"),
])
def test_phase120_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
