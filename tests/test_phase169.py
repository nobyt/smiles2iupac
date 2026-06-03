"""Phase 169: イソシアニド (R-N≡C) の命名

IUPAC 2013 P-62.5.3.2 substitutive PIN:
  [C-]#[N+]C   → isocyanomethane
  [C-]#[N+]CC  → isocyanoethane
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[C-]#[N+]C",   "isocyanomethane"),
    ("[C-]#[N+]CC",  "isocyanoethane"),
    ("[C-]#[N+]CCC", "isocyanopropane"),
    # 回帰: アンモニウムは変わらない
    ("C[NH3+]",         "methylazanium"),
    ("[N+](C)(C)(C)C",  "tetramethylazanium"),
    ("[NH4+]",          "ammonium"),
    # 回帰: ニトリルは変わらない
    ("CC#N", "acetonitrile"),
    ("CCC#N", "propanenitrile"),
])
def test_phase169_isocyanide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
