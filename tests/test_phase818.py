"""Phase 818: imidazo[1,5-b]pyridazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C5 → 5(1H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1ccc2cncn2n1",   "imidazo[1,5-b]pyridazin-2(1H)-one"),
    ("Sc1ccc2cncn2n1",   "imidazo[1,5-b]pyridazin-2(1H)-thione"),
    ("Oc1ncn2ncccc12",   "imidazo[1,5-b]pyridazin-5(1H)-one"),
    ("Sc1ncn2ncccc12",   "imidazo[1,5-b]pyridazin-5(1H)-thione"),
    ("Oc1ncc2cccnn12",   "imidazo[1,5-b]pyridazin-7(1H)-one"),
    ("Sc1ncc2cccnn12",   "imidazo[1,5-b]pyridazin-7(1H)-thione"),
    ("c1cnn2cncc2c1",    "imidazo[1,5-b]pyridazine"),
])
def test_phase818(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
