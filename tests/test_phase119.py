"""Phase 119: 2炭素鎖 amine/thiol のロカント省略 (IUPAC 2013 P-31.1.2.2)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ethanamine (not ethan-1-amine)
    ("CCN", "ethanamine"),
    # ethanethiol (not ethane-1-thiol)
    ("CCS", "ethanethiol"),
    # 回帰: 3炭素以上はロカント保持
    ("CCCN", "propan-1-amine"),
    ("CC(N)C", "propan-2-amine"),
    ("CCCS", "propane-1-thiol"),
    ("CC(S)C", "propane-2-thiol"),
    # 回帰: 1炭素はロカント省略 (既存動作)
    ("CN", "methanamine"),
    ("CS", "methanethiol"),
])
def test_phase119_two_carbon_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
