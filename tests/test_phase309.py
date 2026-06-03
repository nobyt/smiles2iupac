"""Phase 309: benzene dicarbonitrile naming (IUPAC 2013 P-66.5.1.1.1).

Two nitrile groups on benzene → "benzene-X,Y-dicarbonitrile".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("N#Cc1ccccc1C#N",   "benzene-1,2-dicarbonitrile"),
    ("N#Cc1cccc(C#N)c1", "benzene-1,3-dicarbonitrile"),
    ("N#Cc1ccc(C#N)cc1", "benzene-1,4-dicarbonitrile"),
    # regressions: chain dinitriles unchanged
    ("N#CCCC#N",         "butanedinitrile"),
    ("N#CCC#N",          "propanedinitrile"),
    # regression: mono benzonitrile unchanged
    ("N#Cc1ccccc1",      "benzonitrile"),
])
def test_phase309_benzene_dicarbonitrile(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
