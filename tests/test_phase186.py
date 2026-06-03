"""Phase 186: N'-置換ヒドラジド (IUPAC 2013 P-66.3.5)

  CC(=O)NNC  → N'-methylethanohydrazide   (N' = terminal N, 末端側)
  CC(=O)NNCC → N'-ethylethanohydrazide
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N'-置換 (terminal N 側)
    ("CC(=O)NNC",    "N'-methylethanohydrazide"),
    ("CC(=O)NNCC",   "N'-ethylethanohydrazide"),
    ("CCCC(=O)NNC",  "N'-methylbutanohydrazide"),
    # 無置換 (回帰)
    ("CC(=O)NN",     "ethanohydrazide"),
    ("CCCC(=O)NN",   "butanohydrazide"),
    ("CCC(=O)NN",    "propanohydrazide"),
    # 芳香族 (回帰)
    ("O=C(NN)c1ccccc1", "benzohydrazide"),
])
def test_phase186_n_substituted_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
