"""Phase 208: cyanamide retained name (IUPAC 2013 P-66.4.1.1)

  NC#N   → cyanamide
  CNC#N  → N-methylcyanamide
  CN(C)C#N → N,N-dimethylcyanamide

H₂N-C≡N is cyanamide; N-substituted derivatives retain the parent name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Parent
    ("NC#N",       "cyanamide"),
    # N-monosubstituted
    ("CNC#N",      "N-methylcyanamide"),
    ("CCNC#N",     "N-ethylcyanamide"),
    # N,N-disubstituted
    ("CN(C)C#N",   "N,N-dimethylcyanamide"),
    # regression: simple nitrile unaffected
    ("CC#N",       "acetonitrile"),
    ("CCC#N",      "propanenitrile"),
])
def test_phase208_cyanamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
