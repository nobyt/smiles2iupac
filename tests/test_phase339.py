"""Phase 339: E/Z in substituted hydrazone, semicarbazone, and thiosemicarbazone chains (IUPAC 2013).

When the parent carbonyl chain of an N-substituted hydrazone, semicarbazone, or
thiosemicarbazone contains a C=C double bond with defined E/Z geometry, the stereo
descriptor and unsaturation locant must appear in the name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-substituted hydrazone with E/Z chain
    ("C/C=C/CC(=NNC)C",              "(4E)-hex-4-en-2-one N-methylhydrazone"),
    (r"C/C=C\CC(=NNC)C",             "(4Z)-hex-4-en-2-one N-methylhydrazone"),
    ("C/C=C/CC(=NNCC)C",             "(4E)-hex-4-en-2-one N-ethylhydrazone"),
    # aldehyde N-substituted hydrazone with E/Z chain
    ("C/C=C/CC=NNC",                 "(3E)-pent-3-enal N-methylhydrazone"),
    # semicarbazone with E/Z chain
    ("C/C=C/CC(=NNC(=O)N)C",         "(4E)-hex-4-en-2-one semicarbazone"),
    (r"C/C=C\CC(=NNC(=O)N)C",        "(4Z)-hex-4-en-2-one semicarbazone"),
    # thiosemicarbazone with E/Z chain
    ("C/C=C/CC(=NNC(=S)N)C",         "(4E)-hex-4-en-2-one thiosemicarbazone"),
    # regressions: saturated chains unchanged
    ("CC(=NNC)C",                    "propan-2-one N-methylhydrazone"),
    ("CC=NNC",                       "ethanal N-methylhydrazone"),
    ("CC(=NNC(N)=O)C",               "propan-2-one semicarbazone"),
    ("CC(=NNC(N)=S)C",               "propan-2-one thiosemicarbazone"),
    # regression: unsubstituted hydrazone E/Z (was already working)
    ("C/C=C/CC(=NN)C",               "(4E)-hex-4-en-2-one hydrazone"),
])
def test_phase339_ez_substituted_hydrazone_semicarbazone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
