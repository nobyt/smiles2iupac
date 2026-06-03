"""Phase 172: 環状化合物のロカント 1 省略 (IUPAC 2013 P-63.5.1.1)

置換基・多重結合がない場合、ロカント 1 を省略:
  OC1CCCCC1  → cyclohexanol    (not cyclohexan-1-ol)
  O=C1CCCCC1 → cyclohexanone   (not cyclohexan-1-one)
  C1=CCCCC1  → cyclohexene     (not cyclohex-1-ene)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # アルコール: ロカント 1 省略
    ("OC1CCCCC1",  "cyclohexanol"),
    ("OC1CCCC1",   "cyclopentanol"),
    ("OC1CC1",     "cyclopropanol"),
    # ケトン
    ("O=C1CCCCC1", "cyclohexanone"),
    ("O=C1CCCC1",  "cyclopentanone"),
    # チオール
    ("SC1CCCCC1",  "cyclohexanethiol"),
    ("SC1CCCC1",   "cyclopentanethiol"),
    # アルケン
    ("C1=CCCCC1",  "cyclohexene"),
    ("C1=CCCC1",   "cyclopentene"),
    ("C1=CCC1",    "cyclobutene"),
    # 置換体はロカントが必要
    ("OC1C=CCCC1", "cyclohex-2-en-1-ol"),
    ("CC1CCCCC1O", "2-methylcyclohexan-1-ol"),
    # 回帰: 置換基なし環は変わらない
    ("C1CCCCC1",   "cyclohexane"),
    ("c1ccccc1",   "benzene"),
])
def test_phase172_cyclic_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
