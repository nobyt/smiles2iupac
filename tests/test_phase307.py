"""Phase 307: benzene dicarboxamide naming (IUPAC 2013 P-65.1.2.4).

Two amide groups on benzene → "benzene-X,Y-dicarboxamide".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NC(=O)c1ccccc1C(=O)N",   "benzene-1,2-dicarboxamide"),
    ("NC(=O)c1cccc(C(=O)N)c1", "benzene-1,3-dicarboxamide"),
    ("NC(=O)c1ccc(C(=O)N)cc1", "benzene-1,4-dicarboxamide"),
    # regressions: chain diamides unchanged
    ("NC(=O)CC(=O)N",           "propanediamide"),
    ("NC(=O)CCC(=O)N",          "butanediamide"),
    # regression: mono benzamide unchanged
    ("NC(=O)c1ccccc1",          "benzamide"),
])
def test_phase307_benzene_dicarboxamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
