"""Phase 226: polysulfide naming — trisulfide and tetrasulfide (IUPAC 2013).

R-S-S-S-R' → dialkyl trisulfide; R-S-S-S-S-R' → dialkyl tetrasulfide.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # trisulfide
    ("CSSSC",    "dimethyl trisulfide"),
    ("CCSSSC",   "ethyl methyl trisulfide"),
    ("CCSSSCC",  "diethyl trisulfide"),
    # tetrasulfide
    ("CSSSSC",   "dimethyl tetrasulfide"),
    # regression: disulfide still works
    ("CSSC",     "dimethyl disulfide"),
    ("CCSSCC",   "diethyl disulfide"),
])
def test_phase226_polysulfide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
