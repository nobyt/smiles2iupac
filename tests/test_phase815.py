"""Phase 815: pyrazolo[1,5-b][1,2,4]triazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C3 → 3(3H)-one/thione; C7 → 7(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cnn2nccc2n1",   "pyrazolo[1,5-b][1,2,4]triazin-2(1H)-one"),
    ("Sc1cnn2nccc2n1",   "pyrazolo[1,5-b][1,2,4]triazin-2(1H)-thione"),
    ("Oc1cnc2ccnn2n1",   "pyrazolo[1,5-b][1,2,4]triazin-3(3H)-one"),
    ("Sc1cnc2ccnn2n1",   "pyrazolo[1,5-b][1,2,4]triazin-3(3H)-thione"),
    ("Oc1cc2nccnn2n1",   "pyrazolo[1,5-b][1,2,4]triazin-7(5H)-one"),
    ("Sc1cc2nccnn2n1",   "pyrazolo[1,5-b][1,2,4]triazin-7(5H)-thione"),
    ("c1cnn2nccc2n1",    "pyrazolo[1,5-b][1,2,4]triazine"),
])
def test_phase815(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
