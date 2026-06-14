"""Phase 382: Fix (thiocyanato)benzene and related: N≡C must not be treated as amino.

When -S-C≡N is attached to benzene (canonical SMILES N#CSc1ccccc1), the
nitrile N was incorrectly treated as amino. Fixed and names updated to
functional-class form per Phase 521.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phenyl thiocyanate (the broken case — N appears first in canonical SMILES)
    ("N#CSc1ccccc1",          "phenyl thiocyanate"),
    ("c1ccccc1SC#N",          "phenyl thiocyanate"),
    # alkyl thiocyanates still work
    ("CCSC#N",                "ethyl thiocyanate"),
    ("CSC#N",                 "methyl thiocyanate"),
    # aminomethyl still recognized (real NH2)
    ("NCC",                   "ethanamine"),
    ("NCc1ccccc1",            "phenylmethanamine"),
    # 2-aminoethyl still recognized
    ("NCCC(=O)O",             "3-aminopropanoic acid"),
])
def test_phase382_thiocyanato_benzene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
