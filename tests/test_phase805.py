"""Phase 805: [1,2,4]triazolo[1,5-a]pyrazine α-ol/thiol → tautomers.

- C2 → 2(2H)-one/thione; C5 → 5(1H)-one/thione
- C6 → 6(4H)-one/thione; C8 → 8(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,4]triazolo[1,5-a]pyrazine C2-OH/SH
    ("Oc1nc2cnccn2n1",   "[1,2,4]triazolo[1,5-a]pyrazin-2(2H)-one"),
    ("Sc1nc2cnccn2n1",   "[1,2,4]triazolo[1,5-a]pyrazin-2(2H)-thione"),
    # [1,2,4]triazolo[1,5-a]pyrazine C5-OH/SH
    ("Oc1cncc2ncnn12",   "[1,2,4]triazolo[1,5-a]pyrazin-5(1H)-one"),
    ("Sc1cncc2ncnn12",   "[1,2,4]triazolo[1,5-a]pyrazin-5(1H)-thione"),
    # [1,2,4]triazolo[1,5-a]pyrazine C6-OH/SH
    ("Oc1cn2ncnc2cn1",   "[1,2,4]triazolo[1,5-a]pyrazin-6(4H)-one"),
    ("Sc1cn2ncnc2cn1",   "[1,2,4]triazolo[1,5-a]pyrazin-6(4H)-thione"),
    # [1,2,4]triazolo[1,5-a]pyrazine C8-OH/SH
    ("Oc1nccn2ncnc12",   "[1,2,4]triazolo[1,5-a]pyrazin-8(5H)-one"),
    ("Sc1nccn2ncnc12",   "[1,2,4]triazolo[1,5-a]pyrazin-8(5H)-thione"),
    # Regression: parent ring unchanged
    ("c1cn2ncnc2cn1",    "[1,2,4]triazolo[1,5-a]pyrazine"),
])
def test_phase805_triazolo_pyrazine_15a_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
