"""Phase 402: Barbituric acid → pyrimidine-2,4,6(1H,3H,5H)-trione.

IUPAC 2013 P-31.1.7: for pyrimidine-ring triones, use the aromatic parent
name with indicated-H notation rather than the hexahydro prefix.

Change: first dione block gains an elif for base_name == "hexahydropyrimidine"
with len >= 3 C=O, using all ring atoms with H as indicated-H positions.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # barbituric acid
    ("O=C1NC(=O)NC(=O)C1",       "pyrimidine-2,4,6(1H,3H,5H)-trione"),
    # regression: dihydrouracil (2 C=O) keeps hexahydro name
    ("O=C1NC(=O)NCC1",            "hexahydropyrimidine-2,4-dione"),
    # regression: uracil (aromatic dione) unchanged
    ("O=c1cc[nH]c(=O)[nH]1",     "pyrimidine-2,4(1H,3H)-dione"),
    # regression: thymine unchanged
    ("Cc1c[nH]c(=O)[nH]c1=O",    "5-methylpyrimidine-2,4(1H,3H)-dione"),
    # regression: pyrimidine unchanged
    ("c1ccncn1",                   "pyrimidine"),
    # regression: caprolactam unchanged
    ("O=C1CCCCCN1",                "azepan-2-one"),
])
def test_phase402_barbituric_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
