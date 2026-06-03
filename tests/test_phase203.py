"""Phase 203: alkanediamide naming (IUPAC 2013 P-66.6.3.1)

  NC(=O)C(=O)N  → ethanediamide  (systematic; trivial: oxamide)
  NC(=O)CC(=O)N → propanediamide (trivial: malonamide)
  NC(=O)CCC(=O)N→ butanediamide  (trivial: succinamide)

Two amide groups on a chain yield the suffix -diamide.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NC(=O)C(=O)N",   "ethanediamide"),
    ("NC(=O)CC(=O)N",  "propanediamide"),
    ("NC(=O)CCC(=O)N", "butanediamide"),
    ("NC(=O)CCCC(=O)N","pentanediamide"),
    # regression
    ("NC(=O)N",        "urea"),
    ("CC(=O)N",        "acetamide"),
    ("CCC(=O)N",       "propanamide"),
])
def test_phase203_alkanediamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
