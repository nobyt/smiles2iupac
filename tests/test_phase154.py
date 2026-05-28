"""Phase 154: 飽和5員・6員 O,N / S,N 混合環保留名 (IUPAC 2013)

1,3-oxazolidine, 1,2-oxazolidine, 1,3-thiazolidine,
1,3-oxazinane, 1,3-thiazinane とそれらのラクタム誘導体。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 5員 O,N 環
    ("C1CNCO1",            "1,3-oxazolidine"),
    ("O=C1OCCN1",          "1,3-oxazolidin-2-one"),
    ("O=C1CCNO1",          "1,2-oxazolidin-2-one"),
    # 5員 S,N 環
    ("C1CSCN1",            "1,3-thiazolidine"),
    ("O=C1SCCN1",          "1,3-thiazolidin-2-one"),
    # 6員 O,N 環
    ("O=C1OCCCN1",         "1,3-oxazinan-2-one"),
    # 回帰: 環状ラクタム (Phase 25)
    ("O=C1NCCN1",          "imidazolidin-2-one"),
    ("O=C1OCCCO1",         "1,3-dioxan-2-one"),
    # 回帰: 非環状カルバメート (acyclic carbamate unchanged)
    ("CCOC(=O)N",          "ethyl carbamate"),
    # 回帰
    ("CC(=O)O",            "acetic acid"),
])
def test_phase154_oxazolidine_thiazolidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
