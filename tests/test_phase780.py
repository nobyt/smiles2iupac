"""Phase 780: fused 5+6 bicyclic α-ol/thiol → tautomers (IUPAC 2013).

Systems covered:
- 1H-pyrazolo[3,4-b]pyridine: C3 (alpha to N2) → 3(2H)-one; C6 (alpha to N) → 6(5H)-one
- 1,2-benzoxazole: C3 (alpha to N2) → 3(2H)-one
- imidazo[1,5-a]pyrimidine: C2 → 2(3H)-one; C6 → 6(5H)-one; C8 → 8(5H)-one
- [1,2,3]triazolo[1,5-a]pyrimidine: C5 → 5(4H)-one; C7 → 7(4H)-one
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrazolo[3,4-b]pyridine C3 (alpha to N2)
    ("Oc1n[nH]c2ncccc12",       "1H-pyrazolo[3,4-b]pyridin-3(2H)-one"),
    ("Sc1n[nH]c2ncccc12",       "1H-pyrazolo[3,4-b]pyridin-3(2H)-thione"),
    # 1H-pyrazolo[3,4-b]pyridine C6 (alpha to junction N)
    ("Oc1ccc2cn[nH]c2n1",       "1H-pyrazolo[3,4-b]pyridin-6(5H)-one"),
    ("Sc1ccc2cn[nH]c2n1",       "1H-pyrazolo[3,4-b]pyridin-6(5H)-thione"),
    # 1,2-benzoxazole C3 (alpha to N2)
    ("Oc1noc2ccccc12",           "1,2-benzoxazol-3(2H)-one"),
    ("Sc1noc2ccccc12",           "1,2-benzoxazol-3(2H)-thione"),
    # imidazo[1,5-a]pyrimidine C2 (alpha to N3)
    ("Oc1ccn2cncc2n1",           "imidazo[1,5-a]pyrimidin-2(3H)-one"),
    ("Sc1ccn2cncc2n1",           "imidazo[1,5-a]pyrimidin-2(3H)-thione"),
    # imidazo[1,5-a]pyrimidine C6 (alpha to N5 junction)
    ("Oc1ncc2ncccn12",           "imidazo[1,5-a]pyrimidin-6(5H)-one"),
    ("Sc1ncc2ncccn12",           "imidazo[1,5-a]pyrimidin-6(5H)-thione"),
    # imidazo[1,5-a]pyrimidine C8 (alpha to N5 junction)
    ("Oc1ncn2cccnc12",           "imidazo[1,5-a]pyrimidin-8(5H)-one"),
    ("Sc1ncn2cccnc12",           "imidazo[1,5-a]pyrimidin-8(5H)-thione"),
    # [1,2,3]triazolo[1,5-a]pyrimidine C5 (alpha to N4 junction)
    ("Oc1ccn2nncc2n1",           "[1,2,3]triazolo[1,5-a]pyrimidin-5(4H)-one"),
    ("Sc1ccn2nncc2n1",           "[1,2,3]triazolo[1,5-a]pyrimidin-5(4H)-thione"),
    # [1,2,3]triazolo[1,5-a]pyrimidine C7 (alpha to N4 junction)
    ("Oc1ccnc2cnnn12",           "[1,2,3]triazolo[1,5-a]pyrimidin-7(4H)-one"),
    ("Sc1ccnc2cnnn12",           "[1,2,3]triazolo[1,5-a]pyrimidin-7(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cnc2[nH]ncc2c1",        "1H-pyrazolo[3,4-b]pyridine"),
    ("c1ccc2oncc2c1",            "1,2-benzoxazole"),
    ("c1cnc2cncn2c1",            "imidazo[1,5-a]pyrimidine"),
    ("c1cnc2cnnn2c1",            "[1,2,3]triazolo[1,5-a]pyrimidine"),
])
def test_phase780_fused_bicyclic_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
