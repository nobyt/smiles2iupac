"""Phase 341: E/Z parenthesization in guanidine, nitrosamine, isonitrile (IUPAC 2013).

Fixes parenthesization of E/Z-bearing N-substituents in guanidine and nitrosamine
handlers, and adds E/Z + unsaturation detection to isocyanide substitutive naming.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # guanidine N-substituent parenthesization
    ("C/C=C/CNC(=N)N",             "N-((2E)-but-2-en-1-yl)guanidine"),
    (r"C/C=C\CNC(=N)N",            "N-((2Z)-but-2-en-1-yl)guanidine"),
    # nitrosamine with E/Z chain as parent
    ("C/C=C/CN(N=O)C",             "(2E)-N-methyl-N-nitrosobut-2-en-1-amine"),
    (r"C/C=C\CN(N=O)C",            "(2Z)-N-methyl-N-nitrosobut-2-en-1-amine"),
    # isocyanide E/Z chain
    ("C/C=C/C[N+]#[C-]",           "(2E)-1-isocyanobut-2-ene"),
    (r"C/C=C\C[N+]#[C-]",          "(2Z)-1-isocyanobut-2-ene"),
    # regressions: saturated unchanged
    ("CNC(=N)N",                   "N-methylguanidine"),
    ("CN(N=O)C",                   "N-methyl-N-nitrosomethanamine"),
    ("C[N+]#[C-]",                 "isocyanomethane"),
    ("CC[N+]#[C-]",                "isocyanoethane"),
    ("CC(C)[N+]#[C-]",             "2-isocyanopropane"),
])
def test_phase341_ez_guanidine_nitrosamine_isonitrile(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
