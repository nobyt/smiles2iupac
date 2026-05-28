"""Phase 91: diol/triol/thiol + アルケン鎖の末尾 'e' 保持"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # diol + ene ('e' 保持: but-3-ene-1,2-diol)
    ("C=CC(O)CO", "but-3-ene-1,2-diol"),
    ("C=CCC(O)CO", "pent-4-ene-1,2-diol"),
    ("OCC(O)CC=C", "pent-4-ene-1,2-diol"),
    # thiol + ene ('e' 保持: prop-2-ene-1-thiol)
    ("C=CCS", "prop-2-ene-1-thiol"),
    ("SCC=C", "prop-2-ene-1-thiol"),
    # 回帰: 飽和 diol/thiol
    ("OCC(O)C", "propane-1,2-diol"),
    ("SCCC", "propane-1-thiol"),
])
def test_phase91_diol_thiol_with_ene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
