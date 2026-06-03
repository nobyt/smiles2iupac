"""Phase 127: 環状無水物 (oxa-dione) 命名 (IUPAC 2013 P-31.1.3)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxetane-2,4-dione (malonic anhydride, 4-membered ring)
    ("O=C1CC(=O)O1", "oxetane-2,4-dione"),
    # oxolane-2,5-dione (succinic anhydride, 5-membered ring)
    ("O=C1CCC(=O)O1", "oxolane-2,5-dione"),
    # oxane-2,6-dione (glutaric anhydride, 6-membered ring)
    ("O=C1CCCC(=O)O1", "oxane-2,6-dione"),
    # 回帰: acyclic anhydrides use retained names
    ("CC(=O)OC(=O)C", "acetic anhydride"),
    ("CC(=O)OC(=O)CCC", "acetic butanoic anhydride"),
    # 回帰: lactones unchanged
    ("O=C1CCCO1", "oxolan-2-one"),
    ("O=C1CCCCO1", "oxan-2-one"),
    # 回帰: cyclic imide (systematic PIN)
    ("O=C1NC(=O)CC1", "pyrrolidine-2,5-dione"),
])
def test_phase127_cyclic_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
