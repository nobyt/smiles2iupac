"""Phase 68/175: isothiocyanate substitutive PIN (IUPAC 2013 P-65.5.1.2).

R-N=C=S → isothiocyanato{alkane} (PIN; "{alkyl} isothiocyanate" is retained acceptable).
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CN=C=S",    "isothiocyanatomethane"),
    ("CCN=C=S",   "isothiocyanatoethane"),
    ("CCCN=C=S",  "isothiocyanatopropane"),
    ("CC(C)N=C=S","2-isothiocyanatopropane"),
    ("CCCCN=C=S", "isothiocyanatobutane"),
])
def test_phase68_isothiocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
