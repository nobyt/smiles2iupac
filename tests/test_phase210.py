"""Phase 210: N-substituted guanidine naming (IUPAC 2013 P-66.4.1.4)

  CNC(=N)N    → N-methylguanidine
  CNC(=N)NC   → N,N'-dimethylguanidine
  CN(C)C(=N)N → N,N-dimethylguanidine

Guanidine derivatives retain the 'guanidine' parent with N-substituent prefixes.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Parent retained
    ("NC(=N)N",        "guanidine"),
    # N-monosubstituted
    ("CNC(=N)N",       "N-methylguanidine"),
    ("CCNC(=N)N",      "N-ethylguanidine"),
    # N,N'-disubstituted (one on each N)
    ("CNC(=N)NC",      "N,N'-dimethylguanidine"),
    # N,N-disubstituted (both on same N)
    ("CN(C)C(=N)N",    "N,N-dimethylguanidine"),
    # regression: amidine unaffected
    ("CC(=N)N",        "ethanimidamide"),
    ("CC(=N)NC",       "N-methylethanimidamide"),
])
def test_phase210_substituted_guanidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
