"""Phase 259: salicylic acid, trithiocarbonic acid, amino acid esters (IUPAC 2013).

  OC(=O)c1ccccc1O  → salicylic acid          (retained PIN, IUPAC 2013 P-65.1.1)
  SC(=S)S          → trithiocarbonic acid     (retained)
  NCC(=O)OC        → methyl glycinate         (amino acid ester, P-65.1.2.4)

  Note: nicotinic acid / anthranilic acid / mandelic acid are NOT PINs;
  systematic names are used (pyridine-3-carboxylic acid etc.)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # salicylic acid (retained PIN)
    ("OC(=O)c1ccccc1O", "salicylic acid"),
    # trithiocarbonic acid
    ("SC(=S)S",         "trithiocarbonic acid"),
    # amino acid esters (retained names, IUPAC 2013 P-65.1.2.4)
    ("NCC(=O)OC",       "methyl glycinate"),
    ("NCC(=O)OCC",      "ethyl glycinate"),
    ("CC(N)C(=O)OC",    "methyl alaninate"),
    ("CC(N)C(=O)OCC",   "ethyl alaninate"),
    # regression: amino acids unchanged
    ("NCC(=O)O",        "glycine"),
    ("CC(N)C(=O)O",     "alanine"),
    # regression: benzoic acid unchanged
    ("OC(=O)c1ccccc1",  "benzoic acid"),
    # regression: pyridine carboxylic acid uses systematic name (not nicotinic acid)
    ("OC(=O)c1cccnc1",  "pyridine-3-carboxylic acid"),
])
def test_phase259_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
