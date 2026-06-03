"""Phase 192: イミデートエステル命名 (IUPAC 2013 P-65.1.2.4)

  CC(=N)OCC → ethyl ethanimidate
  C(=N)OC   → methyl methanimidate

構造: C(=N)(O-R) — イミド酸 C(=NH)OH のエステル。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # イミデートエステル
    ("C(=N)OC",           "methyl methanimidate"),
    ("CC(=N)OC",          "methyl ethanimidate"),
    ("CC(=N)OCC",         "ethyl ethanimidate"),
    ("CCC(=N)OC",         "methyl propanimidate"),
    ("CCC(=N)OCC",        "ethyl propanimidate"),
    # 回帰: イミンは変わらない
    ("CC=N",              "ethanimine"),
    ("CC(=N)C",           "propan-2-imine"),
    # 回帰: エステルは変わらない
    ("CC(=O)OC",          "methyl acetate"),
    ("CCC(=O)OCC",        "ethyl propanoate"),
    # 回帰: アミドは変わらない
    ("CC(=O)N",           "acetamide"),
])
def test_phase192_imidate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
