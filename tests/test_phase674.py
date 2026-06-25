"""Phase 674: indoline and isoindoline methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indoline (2,3-dihydro-1H-indole): N(1,H)-C(2)-C(3)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 133; N1-H methylable)
    ("c1ccc2c(c1)CCN2",    "indoline"),
    ("CN1CCc2ccccc21",     "1-methylindoline"),
    ("c1ccc2c(c1)CC(C)N2", "2-methylindoline"),
    ("CC1CNc2ccccc21",     "3-methylindoline"),
    ("Cc1cccc2c1CCN2",     "4-methylindoline"),
    ("Cc1ccc2c(c1)CCN2",   "5-methylindoline"),
    ("Cc1ccc2c(c1)NCC2",   "6-methylindoline"),
    ("Cc1cccc2c1NCC2",     "7-methylindoline"),
    # isoindoline (2,3-dihydro-1H-isoindole): C(1)-N(2,H)-C(3)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 133; C2v-symmetric: 1≡3, 4≡7, 5≡6; N2-H methylable)
    ("c1ccc2c(c1)CNC2",    "isoindoline"),
    ("c1ccc2c(c1)C(C)NC2", "1-methylisoindoline"),
    ("c1ccc2c(c1)CN(C)C2", "2-methylisoindoline"),
    ("c1cc(C)c2c(c1)CNC2", "4-methylisoindoline"),
    ("Cc1ccc2c(c1)CNC2",   "5-methylisoindoline"),
])
def test_phase674(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
