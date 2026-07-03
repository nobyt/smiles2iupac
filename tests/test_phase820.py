"""Phase 820: pyrrolo[1,2-b][1,2,4]triazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C3 → 3(3H)-one/thione; C6 → 6(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cnn2cccc2n1",   "pyrrolo[1,2-b][1,2,4]triazin-2(1H)-one"),
    ("Sc1cnn2cccc2n1",   "pyrrolo[1,2-b][1,2,4]triazin-2(1H)-thione"),
    ("Oc1cnc2cccn2n1",   "pyrrolo[1,2-b][1,2,4]triazin-3(3H)-one"),
    ("Sc1cnc2cccn2n1",   "pyrrolo[1,2-b][1,2,4]triazin-3(3H)-thione"),
    ("Oc1ccc2nccnn12",   "pyrrolo[1,2-b][1,2,4]triazin-6(2H)-one"),
    ("Sc1ccc2nccnn12",   "pyrrolo[1,2-b][1,2,4]triazin-6(2H)-thione"),
    ("c1cc2nccnn2c1",    "pyrrolo[1,2-b][1,2,4]triazine"),
])
def test_phase820(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
