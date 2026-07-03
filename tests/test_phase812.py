"""Phase 812: [1,2,4]triazolo[1,5-b][1,2,4]triazine α-ol/thiol → tautomers.

- C2 → 2(2H)-one/thione; C6 → 6(4H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1nc2nccnn2n1",   "[1,2,4]triazolo[1,5-b][1,2,4]triazin-2(2H)-one"),
    ("Sc1nc2nccnn2n1",   "[1,2,4]triazolo[1,5-b][1,2,4]triazin-2(2H)-thione"),
    ("Oc1cnc2ncnn2n1",   "[1,2,4]triazolo[1,5-b][1,2,4]triazin-6(4H)-one"),
    ("Sc1cnc2ncnn2n1",   "[1,2,4]triazolo[1,5-b][1,2,4]triazin-6(4H)-thione"),
    ("Oc1cnn2ncnc2n1",   "[1,2,4]triazolo[1,5-b][1,2,4]triazin-7(1H)-one"),
    ("Sc1cnn2ncnc2n1",   "[1,2,4]triazolo[1,5-b][1,2,4]triazin-7(1H)-thione"),
    ("c1cnn2ncnc2n1",    "[1,2,4]triazolo[1,5-b][1,2,4]triazine"),
])
def test_phase812(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
