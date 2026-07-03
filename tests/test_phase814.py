"""Phase 814: [1,2,4]triazolo[1,5-b]pyridazine α-ol/thiol → tautomers.

- C2 → 2(2H)-one/thione; C6 → 6(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1nc2cccnn2n1",   "[1,2,4]triazolo[1,5-b]pyridazin-2(2H)-one"),
    ("Sc1nc2cccnn2n1",   "[1,2,4]triazolo[1,5-b]pyridazin-2(2H)-thione"),
    ("Oc1ccc2ncnn2n1",   "[1,2,4]triazolo[1,5-b]pyridazin-6(4H)-one"),
    ("Sc1ccc2ncnn2n1",   "[1,2,4]triazolo[1,5-b]pyridazin-6(4H)-thione"),
    ("c1cnn2ncnc2c1",    "[1,2,4]triazolo[1,5-b]pyridazine"),
])
def test_phase814(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
