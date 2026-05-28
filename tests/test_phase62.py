"""Phase 62: oxo 置換基の修正 (=O → "oxo", not "oxy")"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジオンの dione path
    ("CC(=O)CC(=O)C", "pentane-2,4-dione"),
    ("O=C1CCCC(=O)CC1", "cycloheptane-1,4-dione"),
    # 環状ケトン
    ("O=C1CCCCC1", "cyclohexan-1-one"),
    ("O=C1CCCC1", "cyclopentan-1-one"),
    # アルデヒド + ケトン: アルデヒドが principal, ケトンの O は oxo 置換基
    ("O=CCC(=O)C", "3-oxobutanal"),
    ("O=CCCC(=O)C", "4-oxopentanal"),
    ("O=CCCC(=O)CC", "4-oxohexanal"),
    # カルボン酸 + ケトン: カルボン酸が principal, ケトンの O は oxo 置換基
    ("O=C(O)CCC(C)=O", "4-oxopentanoic acid"),
])
def test_phase62_oxo_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
