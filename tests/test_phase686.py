"""Phase 686: isoindole-1,3(2H)-dione (phthalimide) methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoindole-1,3(2H)-dione: C(1,=O)-N(2,H)-C(3,=O)-C3a-C4-C5-C6-C7-C7a
    # (parent retained as 'phthalimide'; C1/C3=O not methylable; C2v: C4≡C7, C5≡C6)
    ("O=C1NC(=O)c2ccccc21",    "isoindole-1,3(2H)-dione"),
    ("CN1C(=O)c2ccccc2C1=O",   "2-methylisoindole-1,3-dione"),
    ("O=C1NC(=O)c2c(C)cccc21", "4-methylisoindole-1,3(2H)-dione"),
    ("O=C1NC(=O)c2cc(C)ccc21", "5-methylisoindole-1,3(2H)-dione"),
])
def test_phase686(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
