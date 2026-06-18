"""Phase 575: Substituted 9H-pyrido[2,3-b]indole (α-carboline) and
9H-pyrido[3,4-b]indole (β-carboline) naming.
α-carboline: N at 1, NH at 9, C positions 2-8.
β-carboline: N at 2, NH at 9, C positions 1,3-8.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 9H-pyrido[2,3-b]indole (α-carboline)
    ("c1ccc2c(c1)[nH]c1ncccc12",            "9H-pyrido[2,3-b]indole"),
    ("Cc1ccc2c(n1)[nH]c1ccccc12",           "2-methyl-9H-pyrido[2,3-b]indole"),
    ("Cc1cnc2[nH]c3ccccc3c2c1",             "3-methyl-9H-pyrido[2,3-b]indole"),
    ("Cc1ccnc2[nH]c3ccccc3c12",             "4-methyl-9H-pyrido[2,3-b]indole"),
    ("Cc1cccc2[nH]c3ncccc3c12",             "5-methyl-9H-pyrido[2,3-b]indole"),
    ("Cc1ccc2[nH]c3ncccc3c2c1",             "6-methyl-9H-pyrido[2,3-b]indole"),
    ("Cc1ccc2c(c1)[nH]c1ncccc12",           "7-methyl-9H-pyrido[2,3-b]indole"),
    ("Cc1cccc2c1[nH]c1ncccc12",             "8-methyl-9H-pyrido[2,3-b]indole"),
    # 9H-pyrido[3,4-b]indole (β-carboline)
    ("c1ccc2c(c1)[nH]c1cnccc12",            "9H-pyrido[3,4-b]indole"),
    ("Cc1nccc2c1[nH]c1ccccc12",             "1-methyl-9H-pyrido[3,4-b]indole"),
    ("Cc1cc2c(cn1)[nH]c1ccccc12",           "3-methyl-9H-pyrido[3,4-b]indole"),
    ("Cc1cncc2[nH]c3ccccc3c12",             "4-methyl-9H-pyrido[3,4-b]indole"),
    ("Cc1cccc2[nH]c3cnccc3c12",             "5-methyl-9H-pyrido[3,4-b]indole"),
    ("Cc1ccc2[nH]c3cnccc3c2c1",             "6-methyl-9H-pyrido[3,4-b]indole"),
    ("Cc1ccc2c(c1)[nH]c1cnccc12",           "7-methyl-9H-pyrido[3,4-b]indole"),
    ("Cc1cccc2c1[nH]c1cnccc12",             "8-methyl-9H-pyrido[3,4-b]indole"),
])
def test_phase575_carbolines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
