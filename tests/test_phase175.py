"""Phase 175: isocyanate/isothiocyanate functional-class PIN (IUPAC 2013 P-65.3.1).

R-N=C=O → "{alkyl} isocyanate" (PIN):
  CN=C=O   → methyl isocyanate
  CCN=C=O  → ethyl isocyanate
  aromatic: phenyl isocyanate
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 脂肪族
    ("CN=C=O",           "methyl isocyanate"),
    ("CCN=C=O",          "ethyl isocyanate"),
    ("CCCN=C=O",         "propyl isocyanate"),
    ("CCCCN=C=O",        "butyl isocyanate"),
    ("CC(C)N=C=O",       "propan-2-yl isocyanate"),
    # 芳香族
    ("c1ccccc1N=C=O",    "phenyl isocyanate"),
    # 回帰: イソチオシアネート
    ("CN=C=S",           "methyl isothiocyanate"),
    ("CCN=C=S",          "ethyl isothiocyanate"),
])
def test_phase175_isocyanate_functional_class(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
