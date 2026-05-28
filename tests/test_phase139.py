"""Phase 139: 長鎖アルカン/酸の命名 (C16–C30)

IUPAC 2013 P-31.1.1: hexadecane through triacontane
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # C16
    ("CCCCCCCCCCCCCCCC",         "hexadecane"),
    ("CCCCCCCCCCCCCCCC(=O)O",    "hexadecanoic acid"),
    # C17
    ("CCCCCCCCCCCCCCCCC",        "heptadecane"),
    # C18
    ("CCCCCCCCCCCCCCCCCC",       "octadecane"),
    ("CCCCCCCCCCCCCCCCCC(=O)O",  "octadecanoic acid"),
    # C19
    ("CCCCCCCCCCCCCCCCCCC",      "nonadecane"),
    # C20 (already worked)
    ("CCCCCCCCCCCCCCCCCCCC",     "icosane"),
    # C21
    ("CCCCCCCCCCCCCCCCCCCCC",    "henicosane"),
    # C22
    ("CCCCCCCCCCCCCCCCCCCCCC",   "docosane"),
    # C30
    ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "triacontane"),
    # 回帰: C1–C15 unchanged
    ("CC", "ethane"),
    ("CCCCCCCCCCC", "undecane"),
    ("CCCCCCCCCCCCCC", "tetradecane"),
    ("CCCCCCCCCCCCCCC", "pentadecane"),
])
def test_phase139_long_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
