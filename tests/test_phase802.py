"""Phase 802: [1,2,4]triazolo[1,5-a]pyrimidine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C5 → 5(4H)-one/thione; C7 → 7(6H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,4]triazolo[1,5-a]pyrimidine C2-OH/SH
    ("Oc1nc2ncccn2n1",   "[1,2,4]triazolo[1,5-a]pyrimidin-2(1H)-one"),
    ("Sc1nc2ncccn2n1",   "[1,2,4]triazolo[1,5-a]pyrimidin-2(1H)-thione"),
    # [1,2,4]triazolo[1,5-a]pyrimidine C5-OH/SH
    ("Oc1ccn2ncnc2n1",   "[1,2,4]triazolo[1,5-a]pyrimidin-5(4H)-one"),
    ("Sc1ccn2ncnc2n1",   "[1,2,4]triazolo[1,5-a]pyrimidin-5(4H)-thione"),
    # [1,2,4]triazolo[1,5-a]pyrimidine C7-OH/SH
    ("Oc1ccnc2ncnn12",   "[1,2,4]triazolo[1,5-a]pyrimidin-7(6H)-one"),
    ("Sc1ccnc2ncnn12",   "[1,2,4]triazolo[1,5-a]pyrimidin-7(6H)-thione"),
    # Regression: parent ring unchanged
    ("c1cnc2ncnn2c1",    "[1,2,4]triazolo[1,5-a]pyrimidine"),
])
def test_phase802_triazolo_pyrimidine_15a_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
