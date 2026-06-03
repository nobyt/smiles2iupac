"""Phase 248: complete 7- and 8-membered Hantzsch-Widman ring names.

Missing entries: thiepane (7-S), phosphepane (7-P), and all 8-membered
heterocycles (azocane, oxocane, thiocane, phosphocane).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7-membered: 6C + 1 heteroatom
    ("C1CCCCCS1",     "thiepane"),
    ("C1CCCCCP1",     "phosphepane"),
    # 8-membered: 7C + 1 heteroatom
    ("C1CCCCCCN1",    "azocane"),
    ("C1CCCCCCO1",    "oxocane"),
    ("C1CCCCCCS1",    "thiocane"),
    ("C1CCCCCCP1",    "phosphocane"),
    # N-substituted 8-membered ring
    ("CN1CCCCCCC1",   "1-methylazocane"),
    # substituent: thiepan-2-yl
    ("OCC1CCCCCS1",   "(thiepan-2-yl)methanol"),
    # regression: 6-membered unchanged
    ("C1CCCCO1",      "oxane"),
    ("C1CCCCS1",      "thiane"),
    ("C1CCCCCO1",     "oxepane"),
    ("C1CCCCCN1",     "azepane"),
])
def test_phase248_7_8_membered_hw(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
