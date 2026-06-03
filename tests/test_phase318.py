"""Phase 318: cycloalkane dicarbaldehyde, dicarbonitrile, dicarboxamide (IUPAC 2013 P-65.1).

Extends Phase 302/307/309 handlers to support non-aromatic (cycloalkane) rings.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # dicarbaldehyde
    ("O=CC1CCC(C=O)CC1",           "cyclohexane-1,4-dicarbaldehyde"),
    ("O=CC1CCCCC1C=O",             "cyclohexane-1,2-dicarbaldehyde"),
    ("O=CC1CCCC1C=O",              "cyclopentane-1,2-dicarbaldehyde"),
    # dicarbonitrile
    ("N#CC1CCC(C#N)CC1",           "cyclohexane-1,4-dicarbonitrile"),
    ("N#CC1CCCC1C#N",              "cyclopentane-1,2-dicarbonitrile"),
    # dicarboxamide
    ("NC(=O)C1CCC(C(=O)N)CC1",    "cyclohexane-1,4-dicarboxamide"),
    # regressions: benzene di-groups unchanged
    ("O=Cc1ccc(C=O)cc1",           "benzene-1,4-dicarbaldehyde"),
    ("N#Cc1ccc(C#N)cc1",           "benzene-1,4-dicarbonitrile"),
    ("NC(=O)c1ccc(C(=O)N)cc1",    "benzene-1,4-dicarboxamide"),
])
def test_phase318_cycloalkane_digroups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
