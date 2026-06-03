"""Phase 165: ヘテロ芳香環 N-オキシド命名

IUPAC 2013 P-62.4 に基づく命名:
  c1cc[n+]([O-])cc1      → pyridine 1-oxide
  Cc1cc[n+]([O-])cc1     → 4-methylpyridine 1-oxide
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc[n+]([O-])cc1",    "pyridine 1-oxide"),
    ("Cc1cc[n+]([O-])cc1",   "4-methylpyridine 1-oxide"),
    ("Cc1c[n+]([O-])ccn1",   "5-methylpyrazine 1-oxide"),
    # 回帰: 非 N-oxide ヘテロ環は変わらない
    ("c1ccncc1",             "pyridine"),
    ("c1ccoc1",              "furan"),
    ("c1ccc(N=O)cc1",        "nitrosobenzene"),
    # 回帰: 三級アミン N-oxide は変わらない (systematic naming)
    ("C[N+](C)(C)[O-]",      "N,N-dimethylmethanamine N-oxide"),
])
def test_phase165_hetero_n_oxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
