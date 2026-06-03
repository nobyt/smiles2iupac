"""Phase 215: diazonium salt naming (IUPAC 2013 P-68.3.2.4)

  C[N+]#N              → methanediazonium
  CC[N+]#N             → ethanediazonium
  [N+](#N)c1ccccc1     → benzenediazonium

Diazonium compounds R-[N+]#N are named as alkan(arenedia)diazonium ions.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Aliphatic diazonium
    ("C[N+]#N",           "methanediazonium"),
    ("CC[N+]#N",          "ethanediazonium"),
    ("CCC[N+]#N",         "propanediazonium"),
    # Arene diazonium
    ("[N+](#N)c1ccccc1",  "benzenediazonium"),
    # regression: diazo compound unaffected
    ("C=[N+]=[N-]",       "diazomethane"),
    # regression: nitrile unaffected
    ("CC#N",              "acetonitrile"),
])
def test_phase215_diazonium(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
