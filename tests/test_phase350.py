"""Phase 350: O-alkyl carbamothioate (thiocarbamate) naming (IUPAC 2013).

Compounds of the form NH2-C(=S)-O-R or NHR-C(=S)-O-R are named as
O-alkyl carbamothioate or O-alkyl N-substituted carbamothioate.
These are esters of carbamothioic acid (H2N-C(=S)-OH).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Primary amine forms (NH2)
    ("COC(=S)N",                  "O-methyl carbamothioate"),
    ("CCOC(=S)N",                 "O-ethyl carbamothioate"),
    ("CCCOC(=S)N",                "O-propyl carbamothioate"),
    # N-substituted forms
    ("CCOC(=S)NC",                "O-ethyl N-methylcarbamothioate"),
    ("CCOC(=S)NCC",               "O-ethyl N-ethylcarbamothioate"),
    ("CCOC(=S)N(C)C",             "O-ethyl N,N-dimethylcarbamothioate"),
    # regressions: thioamide and carbamate unchanged
    ("CC(=S)N",                   "ethanethioamide"),
    ("CC(=S)NC",                  "N-methylethanethioamide"),
    ("CCOC(=O)N",                 "ethyl carbamate"),
    ("CCOC(=O)NC",                "ethyl N-methylcarbamate"),
])
def test_phase350_o_thiocarbamate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
