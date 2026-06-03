"""Phase 349: E/Z stereo descriptors on C=N bonds (IUPAC 2013).

When a C=N double bond (oxime, imine, hydrazone, semicarbazone, etc.)
carries defined E/Z geometry in the SMILES, the descriptor is prepended
as "(E)-" or "(Z)-".  The lone pair on N is lower priority than any
real substituent, so assignments match the RDKit CIP codes.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Aldehyde oximes
    ("C/C=N/O",                   "(E)-ethanal oxime"),
    (r"C/C=N\O",                  "(Z)-ethanal oxime"),
    ("CCCC/C=N/O",                "(E)-pentanal oxime"),
    # Unstereospecified oximes unchanged
    ("CC=NO",                     "ethanal oxime"),
    ("CC(=NO)C",                  "propan-2-one oxime"),
    # Aldehyde hydrazones
    ("C/C=N/N",                   "(E)-ethanal hydrazone"),
    (r"C/C=N\N",                  "(Z)-ethanal hydrazone"),
    # N-substituted hydrazones
    ("C/C=N/NC",                  "(E)-ethanal N-methylhydrazone"),
    (r"C/C=N\NC",                 "(Z)-ethanal N-methylhydrazone"),
    # Aldehyde semicarbazone
    ("C/C=N/NC(=O)N",             "(E)-ethanal semicarbazone"),
    (r"C/C=N\NC(=O)N",            "(Z)-ethanal semicarbazone"),
    # Aldehyde thiosemicarbazone
    ("C/C=N/NC(=S)N",             "(E)-ethanal thiosemicarbazone"),
    # N-substituted imines
    ("C/C=N/C",                   "(E)-N-methylethanimine"),
    (r"C/C=N\C",                  "(Z)-N-methylethanimine"),
    # O-substituted oximes with E/Z
    ("C/C=N/OC",                  "(E)-O-methylethanal oxime"),
    (r"C/C=N\OC",                 "(Z)-O-methylethanal oxime"),
    # Regressions: unstereospecified unchanged
    ("CC=NC",                     "N-methylethanimine"),
    ("CC=NNC",                    "ethanal N-methylhydrazone"),
    ("CC(=NNC(N)=O)C",            "propan-2-one semicarbazone"),
    ("CC=NOC",                    "O-methylethanal oxime"),
])
def test_phase349_ez_cn_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
