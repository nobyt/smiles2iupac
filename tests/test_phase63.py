"""Phase 63: ジアミン / トリアミン命名"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジアミン (diamine) — 直鎖
    ("NCCN", "ethane-1,2-diamine"),
    ("NCCCN", "propane-1,3-diamine"),
    ("NCCCCN", "butane-1,4-diamine"),
    ("NCCCCCN", "pentane-1,5-diamine"),
    ("NCCCCCCN", "hexane-1,6-diamine"),
    # ジアミン — 分岐鎖
    ("NCC(N)C", "propane-1,2-diamine"),
    ("NCC(N)CC", "butane-1,2-diamine"),
    ("NCC(N)CCC", "pentane-1,2-diamine"),
    # トリアミン
    ("NCCCN(CCCN)CCC", "propane-1,3-triamine"),
])
def test_phase63_diamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
