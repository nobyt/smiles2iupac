"""Phase 233: ethenone (ketene) — drop redundant locant-1 for 2-carbon ene-one.

IUPAC 2013 P-66.6.4.1: locants omitted when unambiguous.
H2C=C=O → ethenone (not ethen-1-one).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C=C=O",   "ethenone"),
    # longer ketenes retain locants
    ("CC=C=O",  "prop-1-en-1-one"),
    ("CCC=C=O", "but-1-en-1-one"),
    # regression: normal ketones unchanged
    ("CC(=O)C",   "acetone"),
    ("CC(=O)CCC", "pentan-2-one"),
])
def test_phase233_ethenone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
