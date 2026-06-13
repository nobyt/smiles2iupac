"""Phase 68/175: isothiocyanate functional-class PIN (IUPAC 2013 P-65.5.1.2).

R-N=C=S → "{alkyl} isothiocyanate" (PIN; e.g., CCN=C=S → ethyl isothiocyanate).
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CN=C=S",    "methyl isothiocyanate"),
    ("CCN=C=S",   "ethyl isothiocyanate"),
    ("CCCN=C=S",  "propyl isothiocyanate"),
    ("CC(C)N=C=S","propan-2-yl isothiocyanate"),
    ("CCCCN=C=S", "butyl isothiocyanate"),
])
def test_phase68_isothiocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
