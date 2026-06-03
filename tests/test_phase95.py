"""Phase 95: 酸無水物 + アルケン鎖 (prop-2-enoic anhydride 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 対称エン酸無水物
    ("C=CC(=O)OC(=O)C=C", "prop-2-enoic anhydride"),
    # 非対称エン+飽和
    ("C=CC(=O)OC(=O)CC", "prop-2-enoic propanoic anhydride"),
    ("C=CC(=O)OC(=O)C", "acetic prop-2-enoic anhydride"),
    # 回帰: 飽和
    ("CC(=O)OC(=O)C", "acetic anhydride"),
    ("CC(=O)OC(=O)CCC", "acetic butanoic anhydride"),
])
def test_phase95_anhydride_with_ene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
