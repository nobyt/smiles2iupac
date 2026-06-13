"""Phase 300: diisocyanate / diisothiocyanate naming (IUPAC 2013 P-65.3.1).

Two isocyanate/isothiocyanate groups on a chain → "X,Y-diisocyanato{alkane}".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # diisocyanate
    ("O=C=NCCN=C=O",   "1,2-diisocyanatoethane"),
    ("O=C=NCCCN=C=O",  "1,3-diisocyanatopropane"),
    ("O=C=NCN=C=O",    "1,1-diisocyanatomethane"),
    # diisothiocyanate
    ("S=C=NCCN=C=S",   "1,2-diisothiocyanatoethane"),
    # regressions: mono forms
    ("O=C=NC",         "methyl isocyanate"),
    ("O=C=NCC",        "ethyl isocyanate"),
    ("S=C=NC",         "methyl isothiocyanate"),
])
def test_phase300_diisocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
