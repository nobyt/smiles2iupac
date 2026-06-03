"""Phase 311: cycloalkane dicarboxylic acid naming (IUPAC 2013 P-65.1.1.1).

Two COOH groups on cycloalkane → "cycloalkane-X,Y-dicarboxylic acid".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("OC(=O)C1CCCCC1C(=O)O",   "cyclohexane-1,2-dicarboxylic acid"),
    ("OC(=O)C1CCC(C(=O)O)CC1", "cyclohexane-1,4-dicarboxylic acid"),
    ("OC(=O)C1CCCC1C(=O)O",    "cyclopentane-1,2-dicarboxylic acid"),
    # regressions: chain diacids unchanged
    ("OC(=O)CCC(=O)O",         "succinic acid"),
    ("OC(=O)CCCC(=O)O",        "glutaric acid"),
    # regressions: benzene diacids (retained names unchanged)
    ("OC(=O)c1ccccc1C(=O)O",   "phthalic acid"),
    # regression: mono acid on ring unchanged
    ("OC(=O)C1CCCCC1",         "cyclohexanecarboxylic acid"),
])
def test_phase311_cycloalkane_diacid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
