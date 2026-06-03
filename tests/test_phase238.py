"""Phase 238: ethenamine — drop redundant locants for 2-carbon vinyl amine.

For C=C-NH2 and N-substituted analogs, the locants "-1-en-1-" are both
unambiguous and are elided to give "ethenamine", "N-methylethenamine", etc.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C=CN",    "ethenamine"),
    ("C=CNC",   "N-methylethenamine"),
    ("C=CNCC",  "N-ethylethenamine"),
    # longer vinyl amines retain locants
    ("CC=CN",   "prop-1-en-1-amine"),
    ("CC=CNC",  "N-methylprop-1-en-1-amine"),
    # regression: aliphatic amines unchanged
    ("CCN",     "ethanamine"),
    ("CCCN",    "propan-1-amine"),
])
def test_phase238_ethenamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
