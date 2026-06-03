"""Phase 302: benzene dicarbaldehyde naming (IUPAC 2013 P-66.6.3.1.1).

Two CHO groups on benzene → "benzene-X,Y-dicarbaldehyde".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=Cc1ccc(C=O)cc1",   "benzene-1,4-dicarbaldehyde"),
    ("O=Cc1ccccc1C=O",     "benzene-1,2-dicarbaldehyde"),
    ("O=Cc1cccc(C=O)c1",   "benzene-1,3-dicarbaldehyde"),
    # regressions: chain dials unchanged
    ("O=CCC=O",            "propanedial"),
    ("O=CCCC=O",           "butanedial"),
    ("O=CCCCC=O",          "pentanedial"),
    # regression: mono benzaldehyde unchanged
    ("O=Cc1ccccc1",        "benzaldehyde"),
])
def test_phase302_benzene_dicarbaldehyde(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
