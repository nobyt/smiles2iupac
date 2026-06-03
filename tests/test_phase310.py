"""Phase 310: benzene dithiol naming (IUPAC 2013 P-63.6.1.2).

Two SH groups on benzene → "benzene-X,Y-dithiol".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Sc1ccccc1S",     "benzene-1,2-dithiol"),
    ("Sc1cccc(S)c1",   "benzene-1,3-dithiol"),
    ("Sc1ccc(S)cc1",   "benzene-1,4-dithiol"),
    # regressions: chain dithiols unchanged
    ("SCCS",           "ethane-1,2-dithiol"),
    ("SCCCS",          "propane-1,3-dithiol"),
    # regression: mono thiophenol unchanged
    ("Sc1ccccc1",      "benzenethiol"),
])
def test_phase310_benzene_dithiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
