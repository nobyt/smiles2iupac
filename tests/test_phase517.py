"""Phase 517: benzeneperoxoic acid — aryl peroxy acid naming (IUPAC 2013 P-65.1.4)

Carbonyl C of peroxy acid directly bonded to benzene ring → benzeneperoxoic acid.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # aryl peroxy acid
    ("c1ccccc1C(=O)OO",  "benzeneperoxoic acid"),
    # aliphatic unaffected
    ("CC(=O)OO",         "ethaneperoxoic acid"),
    ("CCC(=O)OO",        "propaneperoxoic acid"),
    # regression: benzoic acid unaffected
    ("c1ccccc1C(=O)O",   "benzoic acid"),
    # regression: hydroperoxide unaffected
    ("COO",              "methyl hydroperoxide"),
    ("c1ccccc1OO",       "phenyl hydroperoxide"),
])
def test_phase517_benzeneperoxoic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
