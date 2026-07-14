"""Phase 357: Vinyl/alkenyl sulfide, sulfoxide, sulfone naming (IUPAC 2013).

When the parent chain of a sulfide/sulfoxide/sulfone contains a C=C double bond,
the parent name uses the alkene suffix (e.g., "ethene"), not "ethenane".
Previously C=CSC gave "(methylsulfanyl)ethenane" instead of "ethenyl methyl sulfide".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Vinyl sulfone
    ("C=CS(=O)(=O)C",              "ethenyl methyl sulfone"),
    ("C=CS(=O)(=O)CC",             "ethenyl ethyl sulfone"),
    # Vinyl sulfoxide
    ("C=CS(=O)C",                  "ethenyl methyl sulfoxide"),
    # Vinyl sulfide
    ("C=CSC",                      "(methylsulfanyl)ethene"),
    ("C=CSCC",                     "(ethylsulfanyl)ethene"),
    # Internal alkene with sulfonyl/sulfanyl
    ("CC=CS(=O)(=O)C",             "methyl (prop-1-en-1-yl) sulfone"),
    ("CC=CSC",                     "1-(methylsulfanyl)prop-1-ene"),
    # Regressions: saturated cases unchanged
    ("CS(=O)(=O)C",                "dimethyl sulfone"),
    ("CS(=O)C",                    "dimethyl sulfoxide"),
    ("CSC",                        "(methylsulfanyl)methane"),
    ("CCS(=O)(=O)C",               "ethyl methyl sulfone"),
    ("CCSC",                       "(methylsulfanyl)ethane"),
])
def test_phase357_vinyl_sulfide_sulfone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
