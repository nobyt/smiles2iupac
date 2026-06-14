"""Phase 382: Fix (thiocyanato)benzene and related: N≡C must not be treated as amino.

When -S-C≡N is attached to benzene (canonical SMILES N#CSc1ccccc1), the
nitrile N was incorrectly treated as amino, giving '(aminomethylsulfanyl)benzene'
instead of '(thiocyanato)benzene'.

Root cause: _name_carbon_substituent amino-detection used n_h >= 0, which passes
for triple-bond N (n_h=0). Fix: add bond-order check to skip N bonded via triple bond.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phenyl thiocyanate (the broken case — N appears first in canonical SMILES)
    ("N#CSc1ccccc1",          "thiocyanatobenzene"),
    ("c1ccccc1SC#N",          "thiocyanatobenzene"),
    # alkyl thiocyanates still work
    ("CCSC#N",                "thiocyanatoethane"),
    ("CSC#N",                 "thiocyanatomethane"),
    # aminomethyl still recognized (real NH2)
    ("NCC",                   "ethanamine"),
    ("NCc1ccccc1",            "phenylmethanamine"),
    # 2-aminoethyl still recognized
    ("NCCC(=O)O",             "3-aminopropanoic acid"),
])
def test_phase382_thiocyanato_benzene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
