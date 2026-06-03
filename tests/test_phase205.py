"""Phase 205: cyclic diene-dione and diene-diol naming (IUPAC 2013 P-31.1.3)

Ring systems with both multiple bonds and multiple functional groups
need the double bond positions included in the base name.

  O=C1C=CC(=O)C=C1 → cyclohexa-2,5-diene-1,4-dione  (p-benzoquinone)
  O=c1ccc(=O)cc1   → cyclohexa-2,5-diene-1,4-dione  (aromatic SMILES)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # p-benzoquinone: cyclohexa-2,5-diene-1,4-dione
    ("O=C1C=CC(=O)C=C1",  "cyclohexa-2,5-diene-1,4-dione"),
    ("O=c1ccc(=O)cc1",    "cyclohexa-2,5-diene-1,4-dione"),
    # regression: saturated ring diones still work
    ("O=C1CCCCC(=O)C1",   "cycloheptane-1,3-dione"),
    ("O=C1CCC(=O)CC1",    "cyclohexane-1,4-dione"),
    # simple ring ketone regression
    ("O=C1CCCCC1",        "cyclohexanone"),
    ("O=C1CCCC1",         "cyclopentanone"),
])
def test_phase205_cyclic_diene_dione(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
