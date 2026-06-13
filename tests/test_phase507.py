"""Phase 507: purine-2,6-dione derivatives (xanthine, caffeine, theophylline, theobromine)
(IUPAC 2013 P-31.1.7 partially saturated fused heterocycles).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=C1NC=NC2=NC(=O)NC12",       "3,7-dihydro-1H-purine-2,6-dione"),
    ("CN1C(=O)C2=C(N=CN2C)N(C)C1=O", "1,3,7-trimethyl-3,7-dihydro-1H-purine-2,6-dione"),
    ("CN1C=NC2=C1C(=O)NC(=O)N2C",   "1,3-dimethyl-3,7-dihydro-1H-purine-2,6-dione"),
    ("CN1C=NC2=C1N=CNC2=O",         "3,7-dimethyl-3,7-dihydro-1H-purine-2,6-dione"),
    ("O=C1NC=NC2=NC=NC12",          "1,7-dihydro-6H-purin-6-one"),
    ("O=C1N=C2NC(=O)NC2C(=O)N1",    "7,9-dihydro-1H-purine-2,6,8(3H)-trione"),
])
def test_phase507_purine_diones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
