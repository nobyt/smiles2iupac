"""Phase 71: カルバミン酸 (NC(=O)O → carbamic acid)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 非置換
    ("NC(=O)O", "carbamic acid"),
    # N-一置換
    ("CNC(=O)O", "N-methylcarbamic acid"),
    ("CCNC(=O)O", "N-ethylcarbamic acid"),
    ("CCCNC(=O)O", "N-propylcarbamic acid"),
    # N-二置換
    ("CN(C)C(=O)O", "N,N-dimethylcarbamic acid"),
])
def test_phase71_carbamic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
