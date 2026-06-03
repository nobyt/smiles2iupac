"""Phase 80: セミカルバゾン / チオセミカルバゾン命名"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # セミカルバゾン (ketone)
    ("CC(=NNC(N)=O)C", "propan-2-one semicarbazone"),
    ("CCC(=NNC(N)=O)CC", "pentan-3-one semicarbazone"),
    ("CC(=NNC(N)=O)CC", "butan-2-one semicarbazone"),
    # セミカルバゾン (aldehyde)
    ("CC=NNC(N)=O", "ethanal semicarbazone"),
    ("CCC=NNC(N)=O", "propanal semicarbazone"),
    # チオセミカルバゾン
    ("CC(=NNC(N)=S)C", "propan-2-one thiosemicarbazone"),
    ("CC=NNC(N)=S", "ethanal thiosemicarbazone"),
    # ヒドラジドが引き続き正常に動作すること
    ("CC(=O)NN", "ethanohydrazide"),
    ("CC=NN", "ethanal hydrazone"),
])
def test_phase80_semicarbazone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
