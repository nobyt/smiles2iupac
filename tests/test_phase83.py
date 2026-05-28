"""Phase 83: N-置換尿素命名 (N-methylurea, N,N'-dimethylurea 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 片方のN置換
    ("CNC(=O)N", "N-methylurea"),
    ("CCNC(=O)N", "N-ethylurea"),
    ("CCCNC(=O)N", "N-propylurea"),
    # 両方のNに異なる置換基
    ("CNC(=O)NC", "N,N'-dimethylurea"),
    ("CCNC(=O)NCC", "N,N'-diethylurea"),
    ("CNC(=O)NCC", "N-ethyl-N'-methylurea"),
    # 同じNに2つの置換基
    ("CN(C)C(=O)N", "N,N-dimethylurea"),
    # 回帰: 無置換尿素
    ("NC(=O)N", "urea"),
])
def test_phase83_substituted_urea(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
