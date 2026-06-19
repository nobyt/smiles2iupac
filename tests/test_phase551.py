"""Phase 551: Substituted xanthen-9-one and thioxanthen-9-one naming.
Xanthen-9-one (9H-xanthen-9-one) has locants 1-4 on one benzo ring and
5-8 on the other; C9 (carbonyl) and bridge O have no substituent locant.
Hydroxy/amino stay as prefixes (PCG lactone already encoded in the name).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # xanthen-9-one unsubstituted
    ("O=c1c2ccccc2oc2ccccc12",        "xanthen-9-one"),
    # substituted xanthen-9-one
    ("O=c1c2cc(O)ccc2oc2ccccc12",     "2-hydroxyxanthen-9-one"),
    ("O=c1c2ccc(O)cc2oc2ccccc12",     "3-hydroxyxanthen-9-one"),
    ("O=c1c2ccccc2oc2cc(O)ccc12",     "3-hydroxyxanthen-9-one"),
    ("O=c1c2cc(C)ccc2oc2ccccc12",     "2-methylxanthen-9-one"),
    # thioxanthen-9-one unsubstituted
    ("O=c1c2ccccc2sc2ccccc12",        "thioxanthen-9-one"),
    # substituted thioxanthen-9-one
    ("O=c1c2cc(O)ccc2sc2ccccc12",     "2-hydroxythioxanthen-9-one"),
])
def test_phase551_substituted_xanthenone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
