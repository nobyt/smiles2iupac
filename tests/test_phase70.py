"""Phase 70: アミジン (C(=NH)(NH2) → {stem}imidamide)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C(=N)N", "methanimidamide"),
    ("CC(=N)N", "ethanimidamide"),
    ("CCC(=N)N", "propanimidamide"),
    ("CCCC(=N)N", "butanimidamide"),
    ("CCCCC(=N)N", "pentanimidamide"),
])
def test_phase70_amidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
