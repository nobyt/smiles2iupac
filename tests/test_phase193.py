"""Phase 193: イミド酸命名 (IUPAC 2013 P-65.1.2.4)

  CC(=N)O → ethanimidic acid
  C(=N)O  → methanimidic acid

構造: C(=N)-OH — カルボン酸の N-アナログ (アミドの互変異性体)。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # イミド酸
    ("C(=N)O",           "methanimidic acid"),
    ("CC(=N)O",          "ethanimidic acid"),
    ("CCC(=N)O",         "propanimidic acid"),
    ("CCCC(=N)O",        "butanimidic acid"),
    # 回帰: イミデートエステルは変わらない
    ("CC(=N)OCC",        "ethyl ethanimidoate"),
    ("C(=N)OC",          "methyl methanimidoate"),
    # 回帰: イミンは変わらない
    ("CC=N",             "ethan-1-imine"),
    ("CC(=N)C",          "propan-2-imine"),
    # 回帰: カルボン酸は変わらない
    ("CC(=O)O",          "acetic acid"),
    ("CC(=O)N",          "acetamide"),
])
def test_phase193_imidic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
