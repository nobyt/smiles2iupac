"""Phase 670: chromone, xanthen-9-one, and thioxanthen-9-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # chromone (4H-chromen-4-one): O(1)-C(2)-C(3)-C(4,=O)-C4a-...-C8a
    # (parent in Phase 411; C4=O and O1 not methylable)
    # note: ring traversal C4→C3→C2→O1 means atom order 2=C3, 3=C2
    ("O=c1ccoc2ccccc12",    "chromone"),
    ("O=c1cc(C)oc2ccccc12", "2-methylchromone"),
    ("O=c1c(C)coc2ccccc12", "3-methylchromone"),
    ("O=c1ccoc2cccc(C)c12", "5-methylchromone"),
    ("O=c1ccoc2ccc(C)cc12", "6-methylchromone"),
    ("O=c1ccoc2cc(C)ccc12", "7-methylchromone"),
    ("O=c1ccoc2c(C)cccc12", "8-methylchromone"),
    # xanthen-9-one: C9(=O) / O(bridge) / two benzo rings (C2v-symmetric: 1≡8, 2≡7, 3≡6, 4≡5)
    # (parent in Phase 413; C9=O and bridge-O not methylable; 4 unique positions)
    ("O=c1c2ccccc2oc2ccccc12",    "xanthen-9-one"),
    ("O=c1c2c(C)cccc2oc2ccccc12", "1-methylxanthen-9-one"),
    ("O=c1c2cc(C)ccc2oc2ccccc12", "2-methylxanthen-9-one"),
    ("O=c1c2ccc(C)cc2oc2ccccc12", "3-methylxanthen-9-one"),
    ("O=c1c2cccc(C)c2oc2ccccc12", "4-methylxanthen-9-one"),
    # thioxanthen-9-one: same topology with S bridge (C2v-symmetric: 1≡8, 2≡7, 3≡6, 4≡5)
    # (parent in Phase 551; C9=O and bridge-S not methylable; 4 unique positions)
    ("O=c1c2ccccc2sc2ccccc12",    "thioxanthen-9-one"),
    ("O=c1c2c(C)cccc2sc2ccccc12", "1-methylthioxanthen-9-one"),
    ("O=c1c2cc(C)ccc2sc2ccccc12", "2-methylthioxanthen-9-one"),
    ("O=c1c2ccc(C)cc2sc2ccccc12", "3-methylthioxanthen-9-one"),
    ("O=c1c2cccc(C)c2sc2ccccc12", "4-methylthioxanthen-9-one"),
])
def test_phase670(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
