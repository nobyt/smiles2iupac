"""Phase 301: diimine naming (IUPAC 2013 P-62.3.1).

Two imine groups on a chain → "ethane-1,2-diimine" etc.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("N=CC=N",             "ethane-1,2-diimine"),
    ("N=CCCCCC=N",         "hexane-1,6-diimine"),
    ("N=CCCCC=N",          "pentane-1,5-diimine"),
    ("N=CCCC=N",           "butane-1,4-diimine"),
    # regressions: mono imine unchanged
    ("CC(=N)C",            "propan-2-imine"),
    ("CC=N",               "ethanimine"),
])
def test_phase301_diimine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
