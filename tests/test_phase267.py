"""Phase 267: 1,3-benzodioxole retained name (IUPAC 2013).

  c1ccc2c(c1)OCO2  → 1,3-benzodioxole  (methylenedioxybenzene skeleton)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3-benzodioxole base compound
    ("c1ccc2c(c1)OCO2",         "1,3-benzodioxole"),
    ("C1Oc2ccccc2O1",           "1,3-benzodioxole"),
    # regression: benzofuran and indene unchanged
    ("c1coc2ccccc12",           "benzofuran"),
    ("C1=Cc2ccccc2C1",          "1H-indene"),
])
def test_phase267_benzodioxole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
