"""Phase 185: N-置換アミジン (IUPAC 2013 P-66.4.1.1)

  CC(=N)NC  → N-methylethanimidamide   (N = amine N, 単結合側)
  CC(=NC)N  → N'-methylethanimidamide  (N' = imine N, 二重結合側)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-置換 (amine N 側, single bond)
    ("CC(=N)NC",    "N-methylethanimidamide"),
    ("CC(=N)NCC",   "N-ethylethanimidamide"),
    # N'-置換 (imine N 側, double bond)
    ("CC(=NC)N",    "N'-methylethanimidamide"),
    ("CC(=NCC)N",   "N'-ethylethanimidamide"),
    # N,N'-二置換
    ("CC(=NC)NC",   "N-methyl-N'-methylethanimidamide"),
    # 無置換 (回帰)
    ("CC(=N)N",     "ethanimidamide"),
    ("CCCC(=N)N",   "butanimidamide"),
    # グアニジン (回帰)
    ("NC(=N)N",     "guanidine"),
])
def test_phase185_n_substituted_amidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
