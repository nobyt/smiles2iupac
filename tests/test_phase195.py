"""Phase 195: ニトロソ化合物のロカント命名 (IUPAC 2013 P-61.7.1)

  CCCN=O        → 1-nitrosopropane
  CCC(N=O)C     → 2-nitrosobutane
  CC(N=O)CC     → 2-nitrosobutane

鎖長 ≥ 3 の場合のみロカントを付与。1–2 炭素ではロカント不要。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ロカントあり (鎖長 ≥ 3)
    ("CCCN=O",         "1-nitrosopropane"),
    ("CCCCN=O",        "1-nitrosobutane"),
    ("CCC(N=O)C",      "2-nitrosobutane"),
    ("CC(N=O)CC",      "2-nitrosobutane"),
    ("CCC(N=O)CC",     "3-nitrosopentane"),
    ("CCCC(N=O)C",     "2-nitrosopentane"),
    # ロカントなし (鎖長 1–2)
    ("CN=O",           "nitrosomethane"),
    ("CCN=O",          "nitrosoethane"),
    # 回帰: 芳香族は変わらない
    ("c1ccc(N=O)cc1",  "nitrosobenzene"),
])
def test_phase195_nitroso_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
