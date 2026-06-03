"""Phase 351: S-alkyl carbamothioate and S-alkyl carbamodithioate (IUPAC 2013).

Compounds of the form NH2-C(=O)-S-R are named S-alkyl carbamothioate.
Compounds of the form NH2-C(=S)-S-R are named S-alkyl carbamodithioate.
N-substituted variants use the N-prefix convention.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # S-alkyl carbamothioate (thioester of carbamic acid, S-ester)
    ("CSC(=O)N",                    "S-methyl carbamothioate"),
    ("CCSC(=O)N",                   "S-ethyl carbamothioate"),
    ("CCCSC(=O)N",                  "S-propyl carbamothioate"),
    # N-substituted
    ("CCSC(=O)NC",                  "S-ethyl N-methylcarbamothioate"),
    ("CCSC(=O)NCC",                 "S-ethyl N-ethylcarbamothioate"),
    ("CCSC(=O)N(C)C",               "S-ethyl N,N-dimethylcarbamothioate"),
    # S-alkyl carbamodithioate (both C=S and S-ester)
    ("CSC(=S)N",                    "S-methyl carbamodithioate"),
    ("CCSC(=S)N",                   "S-ethyl carbamodithioate"),
    ("CCSC(=S)NC",                  "S-ethyl N-methylcarbamodithioate"),
    ("CCSC(=S)N(C)C",               "S-ethyl N,N-dimethylcarbamodithioate"),
    # regressions: amide, thioamide, carbamate, O-thiocarbamate unchanged
    ("CC(=O)N",                     "acetamide"),
    ("CC(=S)N",                     "ethanethioamide"),
    ("CCOC(=O)N",                   "ethyl carbamate"),
    ("CCOC(=S)N",                   "O-ethyl carbamothioate"),
    ("CCOC(=O)NC",                  "ethyl N-methylcarbamate"),
])
def test_phase351_s_carbamothioate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
