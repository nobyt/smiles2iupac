"""Phase 201: alkanedithiol naming (IUPAC 2013 P-63.6.1.1)

  SCCS  → ethane-1,2-dithiol
  SCCCS → propane-1,3-dithiol
  SCSS  → propane-1,2-dithiol (gem vs vicinal)

Two thiol groups on an alkane chain are named with the suffix -dithiol.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-dithiols
    ("SCCS",   "ethane-1,2-dithiol"),
    # 1,3-dithiol
    ("SCCCS",  "propane-1,3-dithiol"),
    # 1,4-dithiol
    ("SCCCCS", "butane-1,4-dithiol"),
    # 1,5-dithiol
    ("SCCCCCS", "pentane-1,5-dithiol"),
    # geminal dithiol (1,1-dithiol: though rare)
    # regression: monothiols still work
    ("CS",     "methanethiol"),
    ("CCS",    "ethanethiol"),
    ("CCCS",   "propane-1-thiol"),
])
def test_phase201_alkanedithiol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
