"""Phase 392: Suppress spurious locant '1-' on 1-carbon acyl-halide chains (IUPAC 2013 P-14.5.2).

When the acyl chain has only one carbon (methanoyl), every substituent is
trivially at position 1 — the locant is uninformative and must be omitted.

Bug: ClC(=O)N gave '1-aminomethanoyl chloride' instead of
     'aminomethanoyl chloride'.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # amino on methanoyl (carbamoyl halides)
    ("ClC(=O)N",        "aminomethanoyl chloride"),
    ("BrC(=O)N",        "aminomethanoyl bromide"),
    ("FC(=O)N",         "aminomethanoyl fluoride"),
    # regression: locant still present on 2+ carbon chains
    ("ClC(=O)CN",       "2-aminoethanoyl chloride"),
    ("ClC(=O)CCN",      "3-aminopropanoyl chloride"),
    # regression: no-substituent retained names unchanged
    ("ClC(=O)C",        "acetyl chloride"),
    ("ClC(=O)CC",       "propanoyl chloride"),
    ("ClC(=O)CCC",      "butanoyl chloride"),
])
def test_phase392_methanoyl_no_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
