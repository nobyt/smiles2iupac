"""Phase 164: N-ニトロソアミン (R₂N-N=O) の命名

IUPAC 2013 に基づく置換命名:
  CN(N=O)C   → N-methyl-N-nitrosomethanamine
  CCN(N=O)CC → N-ethyl-N-nitrosoethanamine

また hydrazine 検出から N-N=O パターンを除外する修正も含む。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-ニトロソアミン (対称)
    ("CN(N=O)C",   "N-methyl-N-nitrosomethanamine"),
    ("CCN(N=O)CC", "N-ethyl-N-nitrosoethanamine"),
    # N-ニトロソアミン (非対称)
    ("CN(N=O)CC",  "N-methyl-N-nitrosoethanamine"),
    # 回帰: ヒドラジン命名は変わらない
    ("NN",              "hydrazine"),
    ("NNC",             "methylhydrazine"),
    ("CNNC",            "1,2-dimethylhydrazine"),
    ("CN(N)C",          "1,1-dimethylhydrazine"),
    ("c1ccc(NN)cc1",    "phenylhydrazine"),
    # 回帰: C-ニトロソは変わらない
    ("CCC(N=O)C",       "2-nitrosobutane"),
    ("c1ccc(N=O)cc1",   "nitrosobenzene"),
])
def test_phase164_nitrosamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
