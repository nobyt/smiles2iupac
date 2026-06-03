"""Phase 342: N-substituted sulfamide naming (IUPAC 2013 P-65.3.2).

N-substituted sulfamide compounds (R-NH-S(=O)2-NH2 and analogues) now produce
correct IUPAC substitutive names with N/N' labels assigned in alphabetical order
of the substituent names.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # mono-N-substituted sulfamide
    ("CNS(=O)(=O)N",               "N-methylsulfamide"),
    ("CCNS(=O)(=O)N",              "N-ethylsulfamide"),
    ("CCCNS(=O)(=O)N",             "N-propylsulfamide"),
    # N,N'-disubstituted sulfamide (same substituent)
    ("CNS(=O)(=O)NC",              "N,N'-dimethylsulfamide"),
    ("CCNS(=O)(=O)NCC",            "N,N'-diethylsulfamide"),
    # N,N'-disubstituted sulfamide (different substituents; alphabetical)
    ("CNS(=O)(=O)NCC",             "N-ethyl-N'-methylsulfamide"),
    # N-alkenyl sulfamide with E/Z
    ("C/C=C/CNS(=O)(=O)N",         "N-[(2E)-but-2-en-1-yl]sulfamide"),
    (r"C/C=C\CNS(=O)(=O)N",        "N-[(2Z)-but-2-en-1-yl]sulfamide"),
    # regression: unsubstituted sulfamide unchanged
    ("NS(=O)(=O)N",                "sulfamide"),
])
def test_phase342_n_substituted_sulfamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
