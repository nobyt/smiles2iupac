"""Phase 348: Dialkyl carbonothioate and carbonodithioate naming (IUPAC 2013).

Compounds of the form (RO)2C=S are named as dialkyl carbonothioate.
Compounds of the form (RS)2C=S are named as dialkyl carbonodithioate.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Dialkyl carbonothioate: (RO)2C=S
    ("COC(=S)OC",                "dimethyl carbonothioate"),
    ("CCOC(=S)OCC",              "diethyl carbonothioate"),
    ("CCCOC(=S)OCCC",            "dipropyl carbonothioate"),
    # Mixed alkyl carbonothioate
    ("COC(=S)OCC",               "ethyl methyl carbonothioate"),
    # Dialkyl trithiocarbonate: (RS)2C=S
    ("CSC(=S)SC",                "dimethyl trithiocarbonate"),
    ("CCSC(=S)SCC",              "diethyl trithiocarbonate"),
    # regressions: regular carbonate and thioester unchanged
    ("COC(=O)OC",                "dimethyl carbonate"),
    ("CCOC(=O)OCC",              "diethyl carbonate"),
    ("CSC(=O)C",                 "S-methyl ethanethioate"),
    ("CC(=S)OC",                 "O-methyl ethanethioate"),
])
def test_phase348_carbonothioate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
