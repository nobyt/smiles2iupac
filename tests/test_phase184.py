"""Phase 184: イソシアニド substitutive PIN (IUPAC 2013 P-62.5.3.2)

  C[N+]#[C-]   → isocyanomethane  (substitutive PIN; "methyl isocyanide" is retained)
  CC[N+]#[C-]  → isocyanoethane
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 置換命名 (substitutive PIN per IUPAC 2013 P-62.5.3.2)
    ("C[N+]#[C-]",     "isocyanomethane"),
    ("CC[N+]#[C-]",    "isocyanoethane"),
    ("[C-]#[N+]CCC",   "isocyanopropane"),
    ("[C-]#[N+]CCCC",  "isocyanobutane"),
    # 別の SMILES 表記でも同じ結果
    ("[C-]#[N+]C",     "isocyanomethane"),
    ("[C-]#[N+]CC",    "isocyanoethane"),
    # 回帰: ニトリルは変わらない
    ("CC#N",           "acetonitrile"),
    ("CCC#N",          "propanenitrile"),
    # 回帰: イソシアネート (PIN: functional-class name)
    ("CN=C=O",         "methyl isocyanate"),
    ("CCN=C=O",        "ethyl isocyanate"),
])
def test_phase184_isocyanide_substitutive(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
